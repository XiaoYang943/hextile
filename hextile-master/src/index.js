// import min from 'lodash/min'
// import max from 'lodash/max'
const fs = require('fs');
function equirectangular (center, width, radius = 6371000) {
  const rad2deg = 180 / Math.PI
  const dy = width / radius * rad2deg
  const dx = width / (Math.cos(center[1] / rad2deg) * radius) * rad2deg
  return {
    forward ([lng, lat]) {
      return [(lng - center[0]) / dx, (lat - center[1]) / dy]
    },
    inverse ([x, y]) {
      return [center[0] + x * dx, center[1] + y * dy]
    }
  }
}

function polar2cartesian (theta) {
  return [
    Math.sin(theta / 180 * Math.PI),
    -Math.cos(theta / 180 * Math.PI)
  ]
}

function dotProduct (v1, v2) {
  return v1.reduce((sum, e, i) => sum + e * v2[i], 0)
}

function linearSolver ([alpha1, beta1], [alpha2, beta2]) {
  const DET = alpha1 * beta2 - alpha2 * beta1
  return function (d1, d2) {
    return [
      (beta2 * d1 - beta1 * d2) / DET,
      (-alpha2 * d1 + alpha1 * d2) / DET
    ]
  }
}

function isInside ([lng, lat], linearRing) {
  let isInside = false
  for (let i = 1; i < linearRing.length; i++) {
    const deltaYplus = linearRing[i][1] - lat
    const deltaYminus = lat - linearRing[i - 1][1]
    if (deltaYplus > 0 && deltaYminus <= 0) continue
    if (deltaYplus < 0 && deltaYminus >= 0) continue
    const deltaX = (deltaYplus * linearRing[i - 1][0] + deltaYminus * linearRing[i][0]) /
      (deltaYplus + deltaYminus) - lng
    if (deltaX <= 0) continue
    isInside = !isInside
  }
  return isInside
}

function bbox2geojson (bbox) {
  return {
    type: 'Polygon',
    coordinates: [
      [[bbox[0], bbox[1]], [bbox[2], bbox[1]], [bbox[2], bbox[3]], [bbox[0], bbox[3]], [bbox[0], bbox[1]]]
    ]
  }
}
function max(arr) {
  if (arr.length === 0) {
    throw new Error('数组不能为空');
  }
  let max = arr[0];
  for (let i = 1; i < arr.length; i++) {
    if (arr[i] > max) {
      max = arr[i];
    }
  }
  return max;
}

function min(arr) {
  if (arr.length === 0) {
    throw new Error('数组不能为空');
  }
  let min = arr[0];
  for (let i = 1; i < arr.length; i++) {
    if (arr[i] < min) {
      min = arr[i];
    }
  }
  return min;
}


/**
 * @param {(Object|Object[]|[number, number, number, number])} geojson|bbox - https://tools.ietf.org/html/rfc7946
 * @param {('square'(正方形cell渔网)|'hexagon')} options.shape - default 'square'
 * @param {number} options.width - in metre, default 1000, min 500, max 500000 渔网cell宽高
 * @param {number} options.tilt - in deg, default 0(渔网cell旋转角度)
 * @param {[number, number]} options.center - [lon, lat] of grid origin
 * @param {Object} options.projection
 * @param {Function} options.projection.forward - map lonlat to grid coordinates
 * @param {Function} options.projection.inverse - map grid coordinates to lonlat
 */
const southWest = [103.582, 1.16]
const northEast = [104.1647, 1.48073]
const bbox = [southWest[0], southWest[1], northEast[0], northEast[1]]


let featureCollection = {
  "type": "FeatureCollection",
  "features":[]
}

featureCollection.features = hextile(bbox)

// 将JSON字符串写入文件
fs.writeFile('C:\\Users\\heyiyang\\Desktop\\data.geojson', JSON.stringify(featureCollection), 'utf8', (err) => {
  if (err) {
    console.error('写入文件时发生错误:', err);
  } else {
    console.log('JSON数据已成功写入文件');
  }
});

function hextile (geojson, options = {}) {
  if (
    geojson instanceof Array &&
    geojson.length === 4 &&
    geojson.every(v => typeof v === 'number')
  ) {
    geojson = bbox2geojson(geojson)
  }

  // normalize geojson input
  let input = []
  function extractPolygons (node) {
    if (typeof node !== 'object') return
    if (node instanceof Array) {
      node.forEach(extractPolygons)
    } else if (node.type === 'Polygon') {
      input.push(node.coordinates)
    } else if (node.type === 'MultiPolygon') {
      input.push(...node.coordinates)
    } else if (node.type === 'Feature') {
      extractPolygons(node.geometry)
    } else if (node.type === 'GeometryCollection') {
      node.geometries.forEach(extractPolygons)
    } else if (node.type === 'FeatureCollection') {
      node.features.forEach(extractPolygons)
    }
  }
  extractPolygons(geojson)

  input = input.map(coordinates => ({
    coordinates,
    bbox: [
      min(coordinates[0].map(point => point[0])),
      min(coordinates[0].map(point => point[1])),
      max(coordinates[0].map(point => point[0])),
      max(coordinates[0].map(point => point[1]))
    ]
  }))

  // global bbox
  const bbox = [
    min(input.map(polygon => polygon.bbox[0])),
    min(input.map(polygon => polygon.bbox[1])),
    max(input.map(polygon => polygon.bbox[2])),
    max(input.map(polygon => polygon.bbox[3]))
  ]

  // set default options
  options.shape = options.shape || 'square'
  options.tilt = options.tilt || 0
  options.width = options.width || 1000
  options.width = Math.max(options.width, 500)
  options.width = Math.min(options.width, 500000)
  options.center = options.center || [(bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2]
  options.projection = options.projection || equirectangular(options.center, options.width)

  const forward = options.projection.forward
  const inverse = options.projection.inverse
  // different grid spacing for hexagon
  const step = options.shape === 'hexagon' ? Math.sqrt(3) / 4 : 1
  const grid = {}
  let output = []

  /*
    given grid orientation (beta) and a set of potential endpoints,
    find min and max grid number (not inclusive of extreme values)
  */
  function dRange (beta, endpoints) {
    const dValues = endpoints.map(ep => dotProduct(beta, ep))
    return {
      min: Math.floor(min(dValues) / step + 1),
      max: Math.ceil(max(dValues) / step - 1)
    }
  }

  const corners = [
    forward([bbox[0], bbox[1]]),
    forward([bbox[2], bbox[3]]),
    forward([bbox[0], bbox[3]]),
    forward([bbox[2], bbox[1]])
  ]

  if (options.shape === 'square') {
    const beta0 = polar2cartesian(options.tilt)
    const beta1 = polar2cartesian(options.tilt + 90)
    const dRange0 = dRange(beta0, corners)
    const dRange1 = dRange(beta1, corners)

    // enumerate through all potential grid cell
    for (let i = dRange0.min - 1; i <= dRange0.max + 1; i++) {
      grid[i] = {}
      for (let j = dRange1.min - 1; j <= dRange1.max + 1; j++) {
        grid[i][j] = {}
      }
    }

    input.forEach(polygon => {
      polygon.coordinates.forEach(linearRing => {
        linearRing = linearRing.map(forward)
        for (let n = 0; n < linearRing.length - 1; n++) {
          /*
            If a line segment of the input polygon cuts a grid line,
            label the cell directly above and below the point
            where the line segment cuts the grid line 'keep'.
            This traces out the grid cells lying on the edge of the input polygon
          */
          const beta = [
            linearRing[n + 1][1] - linearRing[n][1],
            linearRing[n][0] - linearRing[n + 1][0]
          ]
          const d = linearRing[n][0] * linearRing[n + 1][1] -
            linearRing[n][1] * linearRing[n + 1][0]

          const iRange = dRange(beta0, [linearRing[n], linearRing[n + 1]])
          const iIntersection = linearSolver(beta0, beta)
          for (let i = iRange.min; i <= iRange.max; i++) {
            const intersection = iIntersection(i, d)
            const j = Math.floor(dotProduct(beta1, intersection))
            grid[i][j].keep = true
            grid[i - 1][j].keep = true
          }

          const jRange = dRange(beta1, [linearRing[n], linearRing[n + 1]])
          const jIntersection = linearSolver(beta1, beta)
          for (let j = jRange.min; j <= jRange.max; j++) {
            const intersection = jIntersection(j, d)
            const i = Math.floor(dotProduct(beta0, intersection))
            grid[i][j].keep = true
            grid[i][j - 1].keep = true
          }
        }
      })
    })

    // translate grid cells into output polygons
    const getIntersection = linearSolver(beta0, beta1)

    for (let _i in grid) {
      const i = +_i
      for (let _j in grid[i]) {
        const j = +_j
        output.push({
          id: [i, j].join('.').replace(/-/g, 'M'),
          type: 'Feature',
          properties: {
            address: [i, j],
            center: inverse(getIntersection(i + 0.5, j + 0.5))
          },
          geometry: {
            type: 'Polygon',
            coordinates: [[
              inverse(getIntersection(i, j)),
              inverse(getIntersection(i + 1, j)),
              inverse(getIntersection(i + 1, j + 1)),
              inverse(getIntersection(i, j + 1))
            ]]
          },
          keep: grid[i][j].keep
        })
      }
    }
  } else if (options.shape === 'hexagon') {
    const beta0 = polar2cartesian(options.tilt)
    const beta1 = polar2cartesian(options.tilt + 60)
    const beta2 = polar2cartesian(options.tilt + 120)
    const dRange0 = dRange(beta0, corners)
    const dRange1 = dRange(beta1, corners)

    for (let i = dRange0.min - 2; i <= dRange0.max + 2; i++) {
      grid[i] = {}
      for (let j = dRange1.min - 2; j <= dRange1.max + 2; j++) {
        grid[i][j] = {
          [j - i + 1]: {},
          [j - i - 1]: {}
        }
      }
    }

    input.forEach(polygon => {
      polygon.coordinates.forEach(linearRing => {
        linearRing = linearRing.map(forward)
        for (let n = 0; n < linearRing.length - 1; n++) {
          // similar logic as above
          const beta = [
            linearRing[n + 1][1] - linearRing[n][1],
            linearRing[n][0] - linearRing[n + 1][0]
          ]
          const d = linearRing[n][0] * linearRing[n + 1][1] -
            linearRing[n][1] * linearRing[n + 1][0]

          const iRange = dRange(beta0, [linearRing[n], linearRing[n + 1]])
          const iIntersection = linearSolver(beta0, beta)
          for (let i = iRange.min; i <= iRange.max; i++) {
            const intersection = iIntersection(i * step, d)
            const j = Math.floor(dotProduct(beta1, intersection) / step)
            const k = Math.floor(dotProduct(beta2, intersection) / step)
            if (i - j + k === 1 || i - j + k === -1) {
              grid[i][j][k].keep = true
              grid[i][j + 1][k + 1].keep = true
            } else {
              grid[i][j + 1][k].keep = true
              grid[i][j][k + 1].keep = true
            }
          }

          const jRange = dRange(beta1, [linearRing[n], linearRing[n + 1]])
          const jIntersection = linearSolver(beta1, beta)
          for (let j = jRange.min; j <= jRange.max; j++) {
            const intersection = jIntersection(j * step, d)
            const i = Math.floor(dotProduct(beta0, intersection) / step)
            const k = Math.floor(dotProduct(beta2, intersection) / step)
            if (i - j + k === 1 || i - j + k === -1) {
              grid[i][j][k].keep = true
              grid[i + 1][j][k + 1].keep = true
            } else {
              grid[i + 1][j][k].keep = true
              grid[i][j][k + 1].keep = true
            }
          }

          const kRange = dRange(beta2, [linearRing[n], linearRing[n + 1]])
          const kIntersection = linearSolver(beta2, beta)
          for (let k = kRange.min; k <= kRange.max; k++) {
            const intersection = kIntersection(k * step, d)
            const i = Math.floor(dotProduct(beta0, intersection) / step)
            const j = Math.floor(dotProduct(beta1, intersection) / step)
            if (i - j + k === 1 || i - j + k === -1) {
              grid[i][j][k].keep = true
              grid[i + 1][j + 1][k].keep = true
            } else {
              grid[i + 1][j][k].keep = true
              grid[i][j + 1][k].keep = true
            }
          }
        }
      })
    })

    const getIntersection = linearSolver(beta0, beta1)

    for (let _i in grid) {
      const i = +_i
      if (!grid[i - 1] || !grid[i + 1]) continue
      for (let _j in grid[i]) {
        const j = +_j
        if (!grid[i][j - 1] || !grid[i][j + 1]) continue
        if ((i + j) % 3 === 0) {
          // combines six triangular grid cells into one hexagon grid cell
          const k = j - i
          output.push({
            id: [i, j].join('.').replace(/-/g, 'M'),
            type: 'Feature',
            properties: {
              address: [i, j],
              center: inverse(getIntersection(i * step, j * step))
            },
            geometry: {
              type: 'Polygon',
              coordinates: [[
                inverse(getIntersection(i * step, (j + 1) * step)),
                inverse(getIntersection((i - 1) * step, j * step)),
                inverse(getIntersection((i - 1) * step, (j - 1) * step)),
                inverse(getIntersection(i * step, (j - 1) * step)),
                inverse(getIntersection((i + 1) * step, j * step)),
                inverse(getIntersection((i + 1) * step, (j + 1) * step))
              ]]
            },
            keep: [
              [i, j, k + 1],
              [i - 1, j, k],
              [i, j - 1, k],
              [i, j, k - 1],
              [i + 1, j, k],
              [i, j + 1, k]
            ].some(([i, j, k]) => {
              return grid[i][j][k].keep
            })
          })
        }
      }
    }
  }

  output = output.filter(feature => {
    if (feature.keep) return true
    /*
      if center of grid cell lies inside one of the input polygons,
      include that grid cell in final output
    */
    const center = feature.properties.center
    return input.some(polygon => {
      if (center[0] < polygon.bbox[0]) return false
      if (center[1] < polygon.bbox[1]) return false
      if (center[0] > polygon.bbox[2]) return false
      if (center[1] > polygon.bbox[3]) return false
      const [outerRing, ...innerRings] = polygon.coordinates
      return isInside(center, outerRing) &&
        innerRings.every(innerRing => !isInside(center, innerRing))
    })
  })

  output.forEach(feature => {
    delete feature.keep
    // copy first point to complete linearRing
    const linearRing = feature.geometry.coordinates[0]
    linearRing.push(linearRing[0])
  })

  return output
}



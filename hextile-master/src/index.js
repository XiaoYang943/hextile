// import min from 'lodash/min'
// import max from 'lodash/max'
const fs = require('fs');

// normalize geojson input
let input = []
/**
 * 等距投影
 * @param center 整个渔网的中心点的经纬度坐标(单位：度、角度)
 * @param width 整个渔网的宽度|高度(单位：米)
 * @param radius 地球半径(单位：米)
 */
function equirectangular (center, width, radius = 6371000) {
  /**
   * 在等距投影中，纬度方向上的距离变化在地图上大致是恒定的
   * width: 使用整个渔网的宽度来代表高度方向上的距离变化，因为在等距投影中，地图的宽度和高度方向的缩放是相等的（假设地图是正方形）。
   * radius: 地球的半径，以米为单位。
   *
   * (width / radius)：计算整个渔网单位高度所对应的地球表面上的弧长。eg:渔网宽度为1000米，则计算一米对应地球表面的弧长(单位：弧度)
   * (width / radius) * (180 / Math.PI)：弧长(弧度制)转纬度(角度制)，即纬度变化量
   * 即dy=计算整个渔网的单位高度所对应的纬度变化量。这个变化量用于经纬度坐标和局部笛卡尔坐标互相转化
   */
  const dy = (width / radius) * (180 / Math.PI) // 地图宽度对应的纬度变化量（以度为单位）
  /**
   * 弧度 = (PI*角度)/180
   * 角度 = (180*弧度)/PI
   *
   * dx:代表了地图纬度方向上经度变化量（单位：度），即
   * 在等距投影中，经度方向上的距离变化不是恒定的，它取决于当前的纬度。
   * 这是因为地球是一个椭球体，而不是一个完美的球体，导致纬线圈的长度随着纬度的增加而减小。因此，我们需要根据当前的中心纬度来调整 dx 的计算。
   *
   * 三角函数的参数通常以弧度为单位，而不是度，所以要把经纬度转为弧度
   *
   * Math.cos 是余弦函数，它接收一个弧度值作为参数，并返回该角的余弦值。
   * 在这里，它计算了中心纬度对应的余弦值。余弦函数在地理和地图学中很有用，因为它可以反映地球表面的曲率
   * 在赤道附近（纬度接近 0°），余弦值接近 1，表明该处的曲率较小，地图上的经度间隔相对均匀。
   * 而在两极附近，余弦值接近 0，表明该处的曲率较大，地图上的经度间隔会被压缩。
   *
   */
  let cosOfTheCenterLatitude = Math.cos(Math.PI * center[1] / 180)  // 整个渔网的中心点的纬度的余弦值，在地球椭球模型上，这个值用于修正因为纬度变化导致的经度方向上距离的变化。
  const dx = width / (cosOfTheCenterLatitude * radius) * (180 / Math.PI)
  return {
    /**
     * 为什么要计算偏移量
     * 1. 中心化坐标系统：
     * 地图投影通常需要一个参考点（或中心点），所有坐标转换都是基于这个点的。通过从输入点的经度中减去整个渔网中心的经度，实际上是在建立一个以整个渔网中心为原点的局部坐标系统
     * 在这个局部坐标系统中，渔网中心的坐标是 [0, 0]，其他点的坐标都是相对于这个中心点的偏移量
     * 2. 简化计算：
     * 在局部坐标系统中进行计算通常比在全球经纬度坐标系统中计算要简单得多
     * 全球经纬度坐标系统是一个复杂的球面坐标系统，而局部坐标系统则通常是一个简单的平面坐标系统。通过将地理坐标转换为局部坐标，我们可以利用平面几何的简便性来进行距离、方向等计算
     */
    // 经纬度坐标转换为局部笛卡尔坐标
    lonLat2cartesian ([lon, lat]) {
      let offsetX = lon - center[0]
      let offsetY = lat - center[1]
      let x = offsetX / dx
      let y = offsetY / dy
      return [x, y]
    },
    // 局部笛卡尔坐标转换为经纬度坐标
    cartesian2lonLat ([x, y]) {
      let offsetX = x * dx
      let offsetY = y * dy
      let lon = center[0] + offsetX
      let lat = center[1] + offsetY
      return [lon,lat]
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

function isInside ([lon, lat], linearRing) {
  let isInside = false
  for (let i = 1; i < linearRing.length; i++) {
    const deltaYplus = linearRing[i][1] - lat
    const deltaYminus = lat - linearRing[i - 1][1]
    if (deltaYplus > 0 && deltaYminus <= 0) continue
    if (deltaYplus < 0 && deltaYminus >= 0) continue
    const deltaX = (deltaYplus * linearRing[i - 1][0] + deltaYminus * linearRing[i][0]) /
      (deltaYplus + deltaYminus) - lon
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
 * @param {('square'(正方形cell渔网)|'hexagon'(正六边形))} options.shape - default 'square'
 * @param {number} options.width - in metre, default 1000, min 500, max 500000 渔网cell宽高(米)
 * @param {number} options.rotationAngle - in deg, default 0(渔网cell旋转角度)
 * @param {[number, number]} options.center - [lon, lat] of grid origin 整个渔网的中心点
 * @param {Object} options.projection
 * @param {Function} options.projection.lonLat2cartesian - map lonlat to grid coordinates
 * @param {Function} options.projection.cartesian2lonLat - map grid coordinates to lonlat
 */
const southWest = [103.582, 1.16]
const northEast = [104.1647, 1.48073]
const bbox = [southWest[0], southWest[1], northEast[0], northEast[1]]

let data1 = {
  "type": "FeatureCollection",
  "features": [
  {
    "type": "Feature",
    "properties": {},
    "geometry": {
      "type": "Polygon",
      "coordinates": [
        [
          [108.620282, 34.283929],
          [108.620282, 34.37588],
          [108.759048, 34.37588],
          [108.759048, 34.283929],
          [108.620282, 34.283929]
        ]
      ]
    }
  }
]
}

let featureCollection = {
  "type": "FeatureCollection",
  "features":[]
}

// bbox data1
featureCollection.features = [...hextile(data1,{
  // shape:'hexagon'
})]
// featureCollection.features = [...hextile(bbox,{
//   // shape:'hexagon'
// }),...bbox.features]

// 将JSON字符串写入文件
fs.writeFile('C:\\Users\\heyiyang\\Desktop\\data.geojson', JSON.stringify(featureCollection), 'utf8', (err) => {
  if (err) {
    console.error('写入文件时发生错误:', err);
  } else {
    console.log('JSON数据已成功写入文件');
  }
});
// 入参是geojson形式的bbox
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
function hextile (geojson, options = {}) {
  // 入参是二维数组形式的bbox
  if (
    geojson instanceof Array &&
    geojson.length === 4 &&
    geojson.every(v => typeof v === 'number')
  ) {
    geojson = bbox2geojson(geojson)
  }

  extractPolygons(geojson)

  /**
   * [[[[108.620282,34.283929],[108.620282,34.37588],[108.759048,34.37588],[108.759048,34.283929],[108.620282,34.283929]]]]
   * =>
   * [{"coordinates":[[[108.620282,34.283929],[108.620282,34.37588],[108.759048,34.37588],[108.759048,34.283929],[108.620282,34.283929]]],"bbox":[108.620282,34.283929,108.759048,34.37588]}]
   */
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
  /**
   * [108.620282,34.283929,108.759048,34.37588]
   */
  const bbox = [
    min(input.map(polygon => polygon.bbox[0])),
    min(input.map(polygon => polygon.bbox[1])),
    max(input.map(polygon => polygon.bbox[2])),
    max(input.map(polygon => polygon.bbox[3]))
  ]

  // 设置参数
  options.shape = options.shape || 'square'
  options.rotationAngle = options.rotationAngle || 0
  options.width = options.width || 1000
  options.width = Math.max(options.width, 500)
  options.width = Math.min(options.width, 500000)
  options.center = options.center || [(bbox[0] + bbox[2]) / 2, (bbox[1] + bbox[3]) / 2]
  options.projection = options.projection || equirectangular(options.center, options.width)

  const lonLat2cartesian = options.projection.lonLat2cartesian
  const cartesian2lonLat = options.projection.cartesian2lonLat
  // different grid spacing for hexagon
  const step = options.shape === 'hexagon' ? Math.sqrt(2) / 4 : 1
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
    lonLat2cartesian([bbox[0], bbox[1]]),
    lonLat2cartesian([bbox[2], bbox[3]]),
    lonLat2cartesian([bbox[0], bbox[3]]),
    lonLat2cartesian([bbox[2], bbox[1]])
  ]

  if (options.shape === 'square') {
    const beta0 = polar2cartesian(options.rotationAngle)
    const beta1 = polar2cartesian(options.rotationAngle + 90)

    const dRange0 = dRange(beta0, corners)
    const dRange1 = dRange(beta1, corners)

    // 枚举所有潜在的网格单元
    for (let i = dRange0.min - 1; i <= dRange0.max + 1; i++) {
      grid[i] = {}
      for (let j = dRange1.min - 1; j <= dRange1.max + 1; j++) {
        grid[i][j] = {}
      }
    }
    input.forEach(polygon => {
      polygon.coordinates.forEach(linearRing => {
        linearRing = linearRing.map(lonLat2cartesian)
        for (let n = 0; n < linearRing.length - 1; n++) {
          /*
            If a line segment of the input polygon cuts a grid line,
            label the cell directly above and below the point
            where the line segment cuts the grid line 'keep'.
            This traces out the grid cells lying on the edge of the input polygon
            如果输入多边形的线段切割网格线，
            在点的正上方和正下方标记单元格
            其中线段切割网格线'keep'
            这会追踪出位于输入多边形边缘的网格单元
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

    // translate grid cells into output polygons 将网格单元转换为输出多边形
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
            center: cartesian2lonLat(getIntersection(i + 0.5, j + 0.5))
          },
          geometry: {
            type: 'Polygon',
            coordinates: [[
              cartesian2lonLat(getIntersection(i, j)),
              cartesian2lonLat(getIntersection(i + 1, j)),
              cartesian2lonLat(getIntersection(i + 1, j + 1)),
              cartesian2lonLat(getIntersection(i, j + 1))
            ]]
          },
          keep: grid[i][j].keep
        })
      }
    }
  } else if (options.shape === 'hexagon') {
    const beta0 = polar2cartesian(options.rotationAngle)
    const beta1 = polar2cartesian(options.rotationAngle + 60)
    const beta2 = polar2cartesian(options.rotationAngle + 120)
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
        linearRing = linearRing.map(lonLat2cartesian)
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
              center: cartesian2lonLat(getIntersection(i * step, j * step))
            },
            geometry: {
              type: 'Polygon',
              coordinates: [[
                cartesian2lonLat(getIntersection(i * step, (j + 1) * step)),
                cartesian2lonLat(getIntersection((i - 1) * step, j * step)),
                cartesian2lonLat(getIntersection((i - 1) * step, (j - 1) * step)),
                cartesian2lonLat(getIntersection(i * step, (j - 1) * step)),
                cartesian2lonLat(getIntersection((i + 1) * step, j * step)),
                cartesian2lonLat(getIntersection((i + 1) * step, (j + 1) * step))
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



export function polar2cartesian (theta) {
  let radian = theta / 180 * Math.PI; // 角度转弧度
  return [
    Math.sin(radian),
    -Math.cos(radian)
  ]
}

export function isInside ([lon, lat], linearRing) {
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

export function bbox2geojson (bbox) {
  return {
    type: 'Polygon',
    coordinates: [
      [[bbox[0], bbox[1]], [bbox[2], bbox[1]], [bbox[2], bbox[3]], [bbox[0], bbox[3]], [bbox[0], bbox[1]]]
    ]
  }
}

export function  max(arr) {
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

export function min(arr) {
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

export function linearSolver ([alpha1, beta1], [alpha2, beta2]) {
  const DET = alpha1 * beta2 - alpha2 * beta1 // x1y2 - x2y1
  return function (d1, d2) {
    return [
      (beta2 * d1 - beta1 * d2) / DET,
      (-alpha2 * d1 + alpha1 * d2) / DET
    ]
  }
}
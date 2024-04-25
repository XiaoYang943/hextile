import FishNet from './FishNet.js'
import fs from "fs";
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

let fishNet = new FishNet(data1)
// let fishNet = new FishNet([108.620282,34.283929,108.759048,34.37588])
let featureCollection = {
  "type": "FeatureCollection",
  "features":[]
}

featureCollection.features = [...fishNet.getFishNet(),...data1.features]

fs.writeFile('C:\\Users\\heyiyang\\Desktop\\data.geojson', JSON.stringify(featureCollection), 'utf8', (err) => {
  if (err) {
    console.error('写入文件时发生错误:', err);
  } else {
    console.log('JSON数据已成功写入文件');
  }
});


'use strict';

const geoRev = require('geo-reverse')

let lat = 0.0;
let lon = 0.0;
let station = '';

function parseArgv() {

  const argvs = process.argv.splice(2);
  const args = {};
  for (let i = 0; i < argvs.length; i++) {
    let key = '';
    let val = '';
    if (argvs[i].startsWith('-')) {
      key = argvs[i].replace(/^-*/, '');
    }

    if (key === '') {
      continue;
    }

    const idx = key.indexOf('=');
    if (idx >= 0) {
      val = key.substring(idx + 1);
      key = key.substring(0, idx);
    } else {
      val = argvs[i + 1];
      i += 1;
    }

    args[key] = val;
  }

  for (const key in args) {
    switch (key) {
      case 'station':
        station = args[key];
        break;
      case 'lat':
        lat = args[key];
        break;
      case 'lon':
        lon = args[key];
        break;
      default:
        break;
    }
  }

}

parseArgv();

const res = geoRev.country( lat, lon);

const code = res[0].isoAlpha3;

console.log(code + ',' + station+',' + lat + ',' + lon);


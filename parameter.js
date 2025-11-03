const RISE_RULE = { right: 3.4, left: 3.4 };
const TWISTS_RAW = {
  right: {
    'cw->cw':   [26.37, 29.94, 33.51],
    'cw->acw':  [29.16, 33.58, 38.00],
    'acw->acw': [26.37, 29.94, 33.51],
    'acw->cw':  [14.20, 16.34, 18.40],
  },
  left: {
    'cw->cw':   [22.65, 29.53, 36.41], //-36.41, -29.53, -22.65
    'cw->acw':  [19.68, 29.28, 38.88], //-38.88, -29.28, -19.68
    'acw->acw': [22.65, 29.53, 36.41], //-36.41, -29.53, -22.65
    'acw->cw':  [21.95, 29.48, 37.01], //-37.01, -29.48, -21.95
  }
};
const TWISTS_ABS = {
  right: Object.fromEntries(Object.entries(TWISTS_RAW.right).map(([k,v])=>[k, v.map(x=>Math.abs(x))])),
  left:  Object.fromEntries(Object.entries(TWISTS_RAW.left ).map(([k,v])=>[k, v.map(x=>Math.abs(x))])),
};
function pairKey(a,b){ return `${a}->${b}`; }
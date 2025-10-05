const RISE_RULE = { right: 3.4, left: 3.4 };
const TWISTS_RAW = {
  right: {
    'cw->cw':   [26.37, 29.94, 33.51],
    'cw->acw':  [29.16, 33.58, 38.00],
    'acw->acw': [26.37, 29.94, 33.51],
    'acw->cw':  [14.2, 16.34, 18.40],
  },
  left: {
    'cw->cw':   [-36.41, -29.53, -22.65],
    'cw->acw':  [-38.88, -29.28, -19.68],
    'acw->acw': [-36.41, -29.53, -22.65],
    'acw->cw':  [-37.01, -29.48, -21.95],
  }
};
const TWISTS_ABS = {
  right: Object.fromEntries(Object.entries(TWISTS_RAW.right).map(([k,v])=>[k, v.map(x=>Math.abs(x))])),
  left:  Object.fromEntries(Object.entries(TWISTS_RAW.left ).map(([k,v])=>[k, v.map(x=>Math.abs(x))])),
};
function pairKey(a,b){ return `${a}->${b}`; }
const topoCodeMap = {
  'parallel': 'P',
  'antiparallel-basket': 'B',
  'antiparallel-chair':  'C',
  'hybrid':              'H',
};
function naInitial(naType) { return naType === 'RNA' ? 'R' : 'D'; }
function naWord(naType)    { return naType === 'RNA' ? 'RNA' : 'DNA'; }
function handShort(hand)   { return hand === 'left' ? 'L' : 'R'; }
function handWord(hand)    { return hand === 'left' ? 'Left' : 'Right'; }
function prefixCode(naType, hand) { return naInitial(naType, hand) + handShort(hand); }
function scopeForLayer(layerIndex1Based, topology) {
  const isAnti = (topology === 'antiparallel-basket' || topology === 'antiparallel-chair');
  if (!isAnti) return 'All';
  return (layerIndex1Based % 2 === 1) ? 'Odd' : 'Even';
}
const _tplCache = new Map();
async function fetchTemplateForLayer(naType, hand, topology, layerIndex1Based, senseCWorACW) {
  const folder = `${naWord(naType)}-${handWord(hand)}`;
  const pref   = prefixCode(naType, hand);
  const topo   = topoCodeMap[topology];
  const handS  = handShort(hand);
  const scope  = scopeForLayer(layerIndex1Based, topology);
  const sense  = (senseCWorACW || 'CW').toUpperCase();
  const tail   = `${naWord(naType)}-${handS}-${topo}-${scope}-${sense}.pdb`;
  const cacheKey = `T:${folder}:${pref}:${tail}`;
  if (_tplCache.has(cacheKey)) return _tplCache.get(cacheKey);
  const tryNums = Array.from({length:12}, (_,i)=>String(i+1).padStart(2,'0'));
  const tried = [];
  for (const nn of tryNums) {
    const fname = `${pref}${nn}.${tail}`;
    const url = `https://raw.githubusercontent.com/S-Ando-Biophysics/4MRNA-G4/main/Templates/${folder}/${fname}`;
    tried.push(url);
    try {
      const res = await fetch(url, { cache:'no-store' });
      if (!res.ok) continue;
      const text = await res.text();
      const parsed = parsePDB(text);
      if (parsed.atoms?.length) {
        _tplCache.set(cacheKey, parsed);
        return parsed;
      }
    } catch(e) {}
  }
  const statusEl = document.getElementById('status');
  if (statusEl) statusEl.textContent = 'Error';
  throw new Error(`Template not found: ${tail}`);
}
const numLayersEl   = document.getElementById("numLayers");
const layersSection = document.getElementById("layersSection");
const gapsSection   = document.getElementById("gapsSection");
const actionsSection= document.getElementById("actionsSection");
const layersGrid    = document.getElementById("layersGrid");
const gapsGrid      = document.getElementById("gapsGrid");
const statusEl      = document.getElementById("status");
const btnGen        = document.getElementById("generateBtn");
const btnDl         = document.getElementById("downloadBtn");
const globalHandEl  = document.getElementById("globalHand");
const topologyEl    = document.getElementById("topology");
const naTypeEl      = document.getElementById("naType");
btnGen.addEventListener('click', async () => {
  btnGen.disabled = true;
  btnDl.disabled  = true;
  if (statusEl) statusEl.textContent = 'Building...';
  showSpinner();
  try {
    resetViewerNow();
    const pdbText = await buildCurrentPDB();
    draw3D(pdbText);
    if (statusEl) statusEl.textContent = 'Done';
    btnDl.disabled = false; 
  } catch (e) {
    console.error(e);
    if (statusEl) statusEl.textContent = 'Error';
    hideSpinner();
  } finally {
    btnGen.disabled = false;
  }
});
function getLayersArray(){ return Array.from(document.querySelectorAll('.layerTemplate')).map(sel => sel.value); }
function makeGlobalHandsArray(){
  const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value||"1")));
  const gaps = Math.max(0, N-1);
  return new Array(gaps).fill(globalHandEl.value || 'right');
}
function applyRiseAndTwistFromRules(L, H, tripleIdxArr){
  const rises = document.querySelectorAll('.gapRise');
  const twists = document.querySelectorAll('.gapTwist');
  for(let i=0;i<H.length;i++){
    const hand = H[i];
    const key = pairKey(L[i], L[i+1]);
    const tri = TWISTS_ABS[hand][key];
    if(!tri) throw new Error(`Twist rule not found for hand=${hand}, key=${key}`);
    const idx = (Array.isArray(tripleIdxArr) ? tripleIdxArr[i] : 1) || 0;
    const twAbs = tri[idx];
    const signedTwist = (hand === 'right' ? +twAbs : -twAbs);
    if(rises[i])  rises[i].value  = RISE_RULE[hand].toFixed(2);
    if(twists[i]) twists[i].value = signedTwist.toFixed(2);
  }
}
function buildControls() {
  const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value||"1")));
  layersGrid.innerHTML = "";
  gapsGrid.innerHTML = "";
  for (let i=0; i<N; i++) {
    const div = document.createElement("div");
    div.className = "gap-card";
    div.innerHTML = `
      <h3>Layer ${i+1}</h3>
      <select data-layer="${i}" class="layerTemplate">
        <option value="cw">Clockwise</option>
        <option value="acw">Anticlockwise</option>
      </select>`;
    layersGrid.appendChild(div);
  }
  for (let i=0; i<N-1; i++) {
    const div = document.createElement("div");
    div.className = "gap-card";
    div.innerHTML = `
      <h3>Layer ${i+1} → Layer ${i+1+1}</h3>
      <div class="grid-gap3">
        <div>
          <label>Rise (Å)</label>
          <input type="number" step="0.01" value="3.4" class="gapRise" data-gap="${i}" />
        </div>
        <div>
          <label>Twist (°)</label>
          <input type="number" step="0.1" value="27.05" class="gapTwist" data-gap="${i}" />
        </div>
      </div>`;
    gapsGrid.appendChild(div);
  }
  layersSection.style.display = "block";
  gapsSection.style.display   = (N>1) ? "block" : "none";
  actionsSection.style.display= "block";
  applyRiseAndTwistFromRules(getLayersArray(), makeGlobalHandsArray(), null);
}
function refreshDefaultsFromGlobal(){
  const L = getLayersArray();
  const H = makeGlobalHandsArray();
  applyRiseAndTwistFromRules(L, H, null);
}
function onTopologyChange(){}
function resSeqForChain(iLayer, N, chainID, topology){
  const rev = (topology === 'antiparallel-basket' && (chainID === 'C' || chainID === 'D'))
           || (topology === 'antiparallel-chair'  && (chainID === 'B' || chainID === 'D'))
           || (topology === 'hybrid'              && (chainID === 'D'));
  return rev ? (N - iLayer) : (iLayer + 1);
}
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
function getLayersArray(){ 
  return Array.from(document.querySelectorAll('.layerTemplate')).map(sel => sel.value); 
}
function makeGlobalHandsArray(){
  const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value||"1")));
  const gaps = Math.max(0, N-1);
  return new Array(gaps).fill(globalHandEl.value || 'right');
}
const DEFAULT_TRIPLE_IDX_BY_KEY = {
  'cw->cw':   1,
  'cw->acw':  1,
  'acw->acw': 1,
  'acw->cw':  1,
};
const ROTATION_DIR = { right: +1, left: -1 };
function applyRiseAndTwistFromRules(L, H, tripleIdxArr){
  const rises  = document.querySelectorAll('.gapRise');
  const twists = document.querySelectorAll('.gapTwist');
  for (let i = 0; i < H.length; i++) {
    const hand = H[i];
    const key  = pairKey(L[i], L[i+1]);
    let tri = (typeof TWISTS_RAW !== 'undefined' &&
               TWISTS_RAW[hand] &&
               TWISTS_RAW[hand][key])
      ? TWISTS_RAW[hand][key]
      : null;
    const usingABS = !tri &&
                     typeof TWISTS_ABS !== 'undefined' &&
                     TWISTS_ABS[hand] &&
                     TWISTS_ABS[hand][key];
    if (!tri && usingABS) {
      tri = TWISTS_ABS[hand][key];
    }
    if (!tri) {
      throw new Error(`Twist rule not found for hand=${hand}, key=${key}`);
    }
    const idxFromArg = (Array.isArray(tripleIdxArr) ? tripleIdxArr[i] : undefined);
    const idxDefault = (DEFAULT_TRIPLE_IDX_BY_KEY && DEFAULT_TRIPLE_IDX_BY_KEY[key] != null)
      ? DEFAULT_TRIPLE_IDX_BY_KEY[key]
      : 1;
    const idx = (idxFromArg != null ? idxFromArg : idxDefault) || 0;
    let twistVal = tri[idx];
    if (usingABS) {
      const dir = ROTATION_DIR[hand] ?? +1;
      twistVal = dir * twistVal;
    }
    if (rises[i])  rises[i].value  = Number(RISE_RULE[hand]).toFixed(2);
    if (twists[i]) twists[i].value = Number(twistVal).toFixed(2);
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
  layersGrid.addEventListener('change', (ev) => {
    const el = ev.target;
    if (el && el.classList && el.classList.contains('layerTemplate')) {
      applyRiseAndTwistFromRules(getLayersArray(), makeGlobalHandsArray(), null);
    }
  });
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
numLayersEl.addEventListener('change', buildControls);
globalHandEl.addEventListener('change', refreshDefaultsFromGlobal);
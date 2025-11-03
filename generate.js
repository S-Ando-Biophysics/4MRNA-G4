let _CURRENT_NA_TYPE = 'DNA';
function residueEnergy(resAtoms){
  const ideals = IDEALS_DATA[_CURRENT_NA_TYPE] || IDEALS_DATA.DNA;
  let E=0;
  for(const [n1,n2,d0,w] of ideals.BONDS){
    const d = bondLen(resAtoms,n1,n2);
    if(d!=null){ const e = d-d0; E += (w||1)*e*e; }
  }
  for(const [n1,n2,n3,a0,w] of ideals.ANGLES){
    const ang = angleDeg(resAtoms,n1,n2,n3);
    if(ang!=null){ const e = ang-a0; E += (w||0.5)*e*e; }
  }
  return E;
}
function trySmallRotation(resAtoms, kind, deg){
  if(kind==="c5tors"){
    const ax = c5o5AxisOf(resAtoms); if(!ax) return false;
    const subset = subsetByNames(resAtoms, PHOS_GROUP);
    rotateSubsetAroundAxis(subset, ax.origin, ax.axis, deg*Math.PI/180);
    return true;
  }else if(kind==="glyco"){
    const ax = glycosidicAxisOf(resAtoms); if(!ax) return false;
    for(const a of resAtoms){
      if(!isSugarPhosphateAtom(a)) continue;
      const p={x:a.x,y:a.y,z:a.z};
      const pr = rotateAroundAxis(p, ax.origin, ax.axis, deg*Math.PI/180);
      a.x=pr.x; a.y=pr.y; a.z=pr.z;
    }
    return true;
  }
  return false;
}
function cloneCoords(resAtoms){
  return resAtoms.map(a=>({a, x:a.x,y:a.y,z:a.z}));
}
function restoreCoords(cl){
  for(const o of cl){ o.a.x=o.x; o.a.y=o.y; o.a.z=o.z; }
}
function lsqRefineChain(atoms, chain){
  const seqs = residueSeqListForChain(atoms, chain);
  const kinds = ["c5tors","glyco"];
  let deltas = [2.0, 1.0, 0.5];
  for(const d of deltas){
    for(const rs of seqs){
      const resAtoms = atomsOf(atoms, chain, rs);
      let E0 = residueEnergy(resAtoms);
      let bestGain = 0, bestKind=null, bestSign=0, bestSnap=null;
      for(const kind of kinds){
        for(const sign of [+1,-1]){
          const snap = cloneCoords(resAtoms);
          const ok = trySmallRotation(resAtoms, kind, sign*d);
          if(!ok){ restoreCoords(snap); continue; }
          const E1 = residueEnergy(resAtoms);
          const gain = E0 - E1;
          if(gain > bestGain){ bestGain=gain; bestKind=kind; bestSign=sign; bestSnap=snap; }
          restoreCoords(snap);
        }
      }
      if(bestGain>1e-6){
        trySmallRotation(resAtoms, bestKind, bestSign*d);
      }
    }
  }
}
function runGeometryLSQAllChains(atoms){
  const chains = new Set(atoms.map(a=>a.chainID));
  for(const ch of chains) lsqRefineChain(atoms, ch);
}
function runPOFitTwoPassesForChain(atoms, chain){
  const seqs = residueSeqListForChain(atoms, chain);
  if(seqs.length < 2) return;
  for(let pass=0; pass<2; pass++){
    for(let i=1; i<seqs.length; i++){
      const prev = atomsOf(atoms, chain, seqs[i-1]);
      const curr = atomsOf(atoms, chain, seqs[i]);
      const plan1 = bestAngle_C5torsion(prev, curr);
      apply_C5torsion(curr, plan1);
      const plan2 = bestAngle_Glyco(prev, curr);
      apply_Glyco(curr, plan2);
    }
  }
}
function runPOFitAllChains(atoms){
  const chains = new Set(atoms.map(a=>a.chainID));
  for(const ch of chains) runPOFitTwoPassesForChain(atoms, ch);
}

let outText = "";
async function generatePDB() {
  const naType = (naTypeEl.value === 'RNA') ? 'RNA' : 'DNA';
  _CURRENT_NA_TYPE = naType; 
  const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value || "1")));
  const layerTpls = getLayersArray();
  const rises  = Array.from(document.querySelectorAll(".gapRise")).map(inp => parseFloat(inp.value || "0"));
  const twists = Array.from(document.querySelectorAll(".gapTwist")).map(inp => parseFloat(inp.value || "0"));
  const handGlobal = (globalHandEl.value || 'right');
  const topology = topologyEl?.value || 'parallel';
  const perLayerTemplates = await Promise.all(
    Array.from({length:N}, async (_,i)=>{
      const sense = (layerTpls[i] === 'acw') ? 'ACW' : 'CW';
      return fetchTemplateForLayer(naType, handGlobal, topology, i+1, sense);
    })
  );
  let serial = 1, cumZ = 0.0, cumDeg = 0.0;
  const atomsOut = [];
  const header = [];
  header.push(`REMARK   4MRNA-G4`);
  for (let i = 0; i < N; i++) {
    const raw = perLayerTemplates[i].atoms;
    const layerPlaced = [];
    if (i > 0) {
      const rise  = rises[i - 1] || 0;
      const twist = twists[i - 1] || 0;
      cumZ   += -rise;
      cumDeg += (handGlobal === "right" ? -1 : +1) * twist;
    }
    for (const a0 of raw) {
      const p  = { x: a0.x, y: a0.y, z: a0.z + cumZ };
      let newResName = a0.resName;
      if (naType === 'RNA' && a0.resName && a0.resName.trim() === "DG") newResName = "G";
      const chainID = a0.chainID || 'A';
      const resSeqLayer = resSeqForChain(i, N, chainID, topology);
      layerPlaced.push({
        ...a0,
        x: p.x, y: p.y, z: p.z,
        resSeq: resSeqLayer,
        serial: serial++,
        resName: newResName
      });
    }
    if (cumDeg !== 0) {
      const m = meanC1pXYZ(layerPlaced);
      const originZ = { x: 0, y: 0, z: m ? m.z : 0 };
      const zAxis   = { x: 0, y: 0, z: 1 };
      rotateAtomsAroundAxis(layerPlaced, originZ, zAxis, cumDeg * Math.PI/180);
    }
    atomsOut.push(...layerPlaced);
  }
  runPOFitAllChains(atomsOut);
  runGeometryLSQAllChains(atomsOut);
  runPOFitAllChains(atomsOut);  
  runGeometryLSQAllChains(atomsOut);

  const lines = [...header, ...atomsOut.map(formatAtom), "END"];
  return lines.join("\n");
}

document.addEventListener('DOMContentLoaded', async ()=>{
  const clamp = ()=>{ if(parseInt(numLayersEl.value||'4') > 4) numLayersEl.value = '4'; };
  clamp();
  numLayersEl.addEventListener('change', ()=>{ clamp(); buildControls(); });
  globalHandEl.addEventListener('change', refreshDefaultsFromGlobal);
  topologyEl?.addEventListener('change', onTopologyChange);
  naTypeEl?.addEventListener('change', ()=>{});
  buildControls();
  ensureViewer();
  setTimeout(()=>{ try{ viewer.resize(); }catch(_){} }, 0);
});
document.getElementById("generateBtn").addEventListener("click", async () => {
  statusEl.textContent = '';
  btnDl.disabled = true;
  try {
    outText = await generatePDB();
    const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value||"1")));
    draw3D(outText, N);
    btnDl.disabled = false;
  } catch(e) {
    statusEl.textContent = "Error: " + (e?.message || e);
  }
});
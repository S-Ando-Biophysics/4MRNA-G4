let TARGET_PO = 1.60;
const SUGAR_NAMES = new Set([
  "C1'", "C2'", "C3'", "C4'", "O4'", "C5'", "O5'", "O3'", "O2'",
  "C1*", "C2*", "C3*", "C4*", "O4*", "C5*", "O5*", "O3*", "O2*"
]);
const PHOS_NAMES = new Set(["P","OP1","OP2","OP3","O1P","O2P","O3P"]);
const PHOS_GROUP = new Set(["P","OP1","OP2","OP3","O1P","O2P","O3P","O5'","O5*"]);
const PHOS_AROUND_O5P = new Set(["OP1","OP2","OP3","O1P","O2P","O3P"]);
function isSugarPhosphateAtom(a){
  const nm = a.name.trim().toUpperCase();
  return SUGAR_NAMES.has(nm) || PHOS_NAMES.has(nm);
}
function distance(a,b){ return Math.hypot(a.x-b.x, a.y-b.y, a.z-b.z); }
function atomKey(a){ return `${a.chainID}|${a.resSeq}`; }
function groupByChainRes(atoms){
  const map = new Map();
  for(const a of atoms){
    const key = atomKey(a);
    if(!map.has(key)) map.set(key, []);
    map.get(key).push(a);
  }
  return map;
}
function glycosidicAxisOf(resAtoms){
  const C1p = _pick(resAtoms, "C1'");
  if(!C1p) return null;
  const N9 = _pick(resAtoms, "N9");
  const N1 = _pick(resAtoms, "N1");
  const axisAtom = N9 || N1;
  if(!axisAtom) return null;
  const origin = {x:C1p.x,y:C1p.y,z:C1p.z};
  const axis   = vnorm(vsub({x:axisAtom.x,y:axisAtom.y,z:axisAtom.z}, origin));
  return {origin, axis};
}
function c5o5AxisOf(resAtoms){
  const C5p = _pick(resAtoms, ["C5'","C5*"]);
  const O5p = _pick(resAtoms, ["O5'","O5*"]);
  if(!C5p || !O5p) return null;
  const origin = {x:C5p.x,y:C5p.y,z:C5p.z};
  const axis   = vnorm(vsub({x:O5p.x,y:O5p.y,z:O5p.z}, origin)); 
  return {origin, axis};
}
function o5pAxisOf(resAtoms){
  const O5p = _pick(resAtoms, ["O5'","O5*"]);
  const P   = _pick(resAtoms, ["P"]);
  if(!O5p || !P) return null;
  const origin = {x:O5p.x,y:O5p.y,z:O5p.z};
  const axis   = vnorm(vsub({x:P.x,y:P.y,z:P.z}, origin)); 
  return {origin, axis};
}
function subsetByNames(resAtoms, nameSet){
  const out=[];
  for(const a of resAtoms){
    if(nameSet.has(a.name.trim().toUpperCase())) out.push(a);
  }
  return out;
}
function bestAngle_C5torsion(prevAtoms, currAtoms){
  const O3p = _pick(prevAtoms, ["O3'","O3*"]);
  const P   = _pick(currAtoms, ["P"]);
  const ax  = c5o5AxisOf(currAtoms);
  if(!O3p || !P || !ax) return null;
  let bestDeg = 0, bestCost = Infinity;
  for(let deg=0; deg<360; deg+=3){
    const rad = deg*Math.PI/180;
    const Pp  = rotateAroundAxis({x:P.x,y:P.y,z:P.z}, ax.origin, ax.axis, rad);
    const d = distance(Pp, O3p);
    const cost = (d - TARGET_PO) * (d - TARGET_PO);
    if(cost < bestCost){ bestCost = cost; bestDeg = deg; }
  }
  let fine = bestDeg, fineCost = bestCost;
  for(let deg=bestDeg-3; deg<=bestDeg+3; deg+=0.2){
    const nd = (deg<0?deg+360:deg%360);
    const rad = nd*Math.PI/180;
    const Pp  = rotateAroundAxis({x:P.x,y:P.y,z:P.z}, ax.origin, ax.axis, rad);
    const d = distance(Pp, O3p);
    const cost = (d - TARGET_PO) * (d - TARGET_PO);
    if(cost < fineCost){ fineCost = cost; fine = nd; }
  }
  return {deg:fine, axis:ax};
}
function apply_C5torsion(currAtoms, plan){
  if(!plan) return;
  const subset = subsetByNames(currAtoms, PHOS_GROUP);
  rotateSubsetAroundAxis(subset, plan.axis.origin, plan.axis.axis, plan.deg*Math.PI/180);
}
function bestAngle_O5torsion(prevAtoms, currAtoms){
  const ax = o5pAxisOf(currAtoms); 
  const O3p = _pick(prevAtoms, ["O3'","O3*"]);
  const P   = _pick(currAtoms, ["P"]);
  if(!ax || !O3p || !P) return null;
  const OP1 = _pick(currAtoms, ["OP1","O1P"]);
  const OP2 = _pick(currAtoms, ["OP2","O2P"]);
  if(!OP1 || !OP2) return null;
  function angleAtP(X){
    const Pp = {x:P.x,y:P.y,z:P.z};
    const v1 = vnorm(vsub(_pick(currAtoms,["O5'","O5*"]), Pp));
    const v2 = vnorm(vsub({x:X.x,y:X.y,z:X.z}, Pp));
    const c = Math.min(1, Math.max(-1, vdot(v1,v2)));
    return Math.acos(c)*180/Math.PI;
  }
  let bestDeg = 0, bestCost = Infinity;
  const subset = subsetByNames(currAtoms, PHOS_AROUND_O5P);
  const snap = subset.map(a=>({a, x:a.x,y:a.y,z:a.z}));
  for(let deg=0; deg<360; deg+=6){
    const rad = deg*Math.PI/180;
    rotateSubsetAroundAxis(subset, ax.origin, ax.axis, rad);
    const a1 = angleAtP(OP1), a2 = angleAtP(OP2);
    const cost = (a1-120)**2 + (a2-120)**2;
    if(cost < bestCost){ bestCost=cost; bestDeg=deg; }
    for(const s of snap){ s.a.x=s.x; s.a.y=s.y; s.a.z=s.z; }
  }
  let fine = bestDeg, fineCost = bestCost;
  for(let deg=bestDeg-6; deg<=bestDeg+6; deg+=0.2){
    const nd = (deg<0?deg+360:deg%360);
    const rad = nd*Math.PI/180;
    rotateSubsetAroundAxis(subset, ax.origin, ax.axis, rad);
    const a1 = angleAtP(OP1), a2 = angleAtP(OP2);
    const cost = (a1-120)**2 + (a2-120)**2;
    if(cost < fineCost){ fineCost=cost; fine=nd; }
    for(const s of snap){ s.a.x=s.x; s.a.y=s.y; s.a.z=s.z; }
  }
  return {deg:fine, axis:ax};
}
function apply_O5torsion(currAtoms, plan){
  if(!plan) return;
  const subset = subsetByNames(currAtoms, PHOS_AROUND_O5P);
  rotateSubsetAroundAxis(subset, plan.axis.origin, plan.axis.axis, plan.deg*Math.PI/180);
}
function bestAngle_Glyco(prevAtoms, currAtoms){
  const O3p = _pick(prevAtoms, ["O3'","O3*"]);
  const P   = _pick(currAtoms, ["P"]);
  const ax  = glycosidicAxisOf(currAtoms);
  if(!O3p || !P || !ax) return null;
  let bestDeg = 0, bestCost = Infinity;
  for(let deg=0; deg<360; deg+=3){
    const rad = deg*Math.PI/180;
    const Pp  = rotateAroundAxis({x:P.x,y:P.y,z:P.z}, ax.origin, ax.axis, rad);
    const d   = distance(Pp, O3p);
    const cost = (d - TARGET_PO) * (d - TARGET_PO);
    if(cost < bestCost){ bestCost = cost; bestDeg = deg; }
  }
  let fine = bestDeg, fineCost = bestCost;
  for(let deg=bestDeg-3; deg<=bestDeg+3; deg+=0.2){
    const nd  = (deg<0?deg+360:deg%360);
    const rad = nd*Math.PI/180;
    const Pp  = rotateAroundAxis({x:P.x,y:P.y,z:P.z}, ax.origin, ax.axis, rad);
    const d   = distance(Pp, O3p);
    const cost = (d - TARGET_PO) * (d - TARGET_PO);
    if(cost < fineCost){ fineCost = cost; fine = nd; }
  }
  return {deg:fine, axis:ax};
}
function apply_Glyco(currAtoms, plan){
  if(!plan) return;
  const origin = plan.axis.origin, axis = plan.axis.axis, rad = plan.deg*Math.PI/180;
  for(const a of currAtoms){
    if(!isSugarPhosphateAtom(a)) continue;
    const p = {x:a.x,y:a.y,z:a.z};
    const pr = rotateAroundAxis(p, origin, axis, rad);
    a.x = pr.x; a.y = pr.y; a.z = pr.z;
  }
}
function residueSeqListForChain(atoms, chain){
  const set = new Set();
  for(const a of atoms){ if(a.chainID===chain) set.add(a.resSeq); }
  return Array.from(set).sort((a,b)=>a-b);
}
function atomsOf(atoms, chain, resSeq){
  const key = `${chain}|${resSeq}`;
  if(!_gbcrCache.has(atoms)) _gbcrCache.set(atoms, groupByChainRes(atoms));
  const byKey = _gbcrCache.get(atoms);
  return (byKey.get(key) || []);
}
const _gbcrCache = new WeakMap();
function detectChainDirection(atoms, chain){
  const seqs = residueSeqListForChain(atoms, chain);
  if (seqs.length < 2) return +1;
  let sumF = 0, sumB = 0, cnt = 0;
  for (let i=1; i<seqs.length; i++){
    const aPrev = atomsOf(atoms, chain, seqs[i-1]);
    const aCurr = atomsOf(atoms, chain, seqs[i]);
    const O3_prev = _pick(aPrev, ["O3'","O3*"]);
    const P_curr  = _pick(aCurr, ["P"]);
    const O3_curr = _pick(aCurr, ["O3'","O3*"]);
    const P_prev  = _pick(aPrev, ["P"]);
    if (O3_prev && P_curr){ sumF += distance(O3_prev, P_curr); cnt++; }
    if (O3_curr && P_prev){ sumB += distance(O3_curr, P_prev); }
  }
  if (cnt === 0) return +1;
  return (sumB < sumF ? -1 : +1);
}
function buildPairsForChain(atoms, chain){
  const seqs = residueSeqListForChain(atoms, chain);
  const dir = detectChainDirection(atoms, chain);
  const pairs = [];
  for (let i=1; i<seqs.length; i++){
    const prev = (dir > 0) ? seqs[i-1] : seqs[i];
    const curr = (dir > 0) ? seqs[i]   : seqs[i-1];
    pairs.push([prev, curr]);
  }
  return { dir, pairs };
}
const IDEALS_DATA = {
  DNA: {
    BONDS: [
      ["C1'", "C2'", 1.52, 1.0],
      ["C1'", "N9",  1.46, 1.0],
      ["C2'", "C3'", 1.52, 1.0],
      ["C3'", "O3'", 1.43, 1.0],
      ["C4'", "C3'", 1.53, 1.0],
      ["C4'", "O4'", 1.45, 1.0],
      ["C5'", "C4'", 1.51, 1.0],
      ["O4'", "C1'", 1.42, 1.0],
      [["O5'","O5*"], "C5'", 1.42, 1.0],
      ["P", ["O5'","O5*"], 1.59, 1.0],
      ["P", "OP1", 1.48, 1.0],
      ["P", "OP2", 1.48, 1.0],
    ],
    ANGLES: [
      ["C1'", "C2'", "C3'", 101.36, 0.5],
      ["C1'", "O4'", "C4'", 109.26, 0.5],
      ["C2'", "C1'", "N9",  114.85, 0.5],
      ["C2'", "C1'", "O4'", 105.54, 0.5],
      ["C2'", "C3'", "C4'", 103.11, 0.5],
      ["C2'", "C3'", "O3'", 111.38, 0.5],
      ["C3'", "C4'", "C5'", 115.08, 0.5],
      ["C3'", "C4'", "O4'", 105.81, 0.5],
      ["C4'", "C3'", "O3'", 108.76, 0.5],
      ["C4'", "C5'", ["O5'","O5*"], 110.04, 0.5],
      ["C5'", "C4'", "O4'", 109.60, 0.5],
      ["C5'", ["O5'","O5*"], "P", 120.29, 0.5],
      [["O5'","O5*"], "P", "OP1", 108.17, 0.5],
      [["O5'","O5*"], "P", "OP2", 107.84, 0.5],
      ["OP1", "P", "OP2", 120.14, 0.5],
      ["N9", "C1'", "O4'", 108.07, 0.5]
    ],
  },
  RNA: {
    BONDS: [
      ["C1'", "C2'", 1.53, 1.0],
      ["C1'", "N9",  1.47, 1.0],
      ["C2'", "C3'", 1.52, 1.0],
      ["C2'", ["O2'","O2*"], 1.42, 1.0],
      ["C3'", "O3'", 1.42, 1.0],
      ["C4'", "C3'", 1.52, 1.0],
      ["C4'", "O4'", 1.45, 1.0],
      ["C5'", "C4'", 1.51, 1.0],
      ["O4'", "C1'", 1.41, 1.0],
      [["O5'","O5*"], "C5'", 1.42, 1.0],
      ["P", ["O5'","O5*"], 1.59, 1.0],
      ["P", "OP1", 1.49, 1.0],
      ["P", "OP2", 1.49, 1.0],
    ],
    ANGLES: [
      ["C1'", "C2'", "C3'", 101.19, 0.5],
      ["C1'", "C2'", ["O2'","O2*"], 109.99, 0.5],
      ["C1'", "O4'", "C4'", 109.54, 0.5],
      ["C2'", "C1'", "N9",  113.68, 0.5],
      ["C2'", "C1'", "O4'", 106.66, 0.5],
      ["C2'", "C3'", "C4'", 102.39, 0.5],
      ["C2'", "C3'", "O3'", 111.62, 0.5],
      ["C3'", "C2'", ["O2'","O2*"], 112.47, 0.5],
      ["C3'", "C4'", "C5'", 115.55, 0.5],
      ["C3'", "C4'", "O4'", 104.73, 0.5],
      ["C4'", "C3'", "O3'", 110.52, 0.5],
      ["C4'", "C5'", ["O5'","O5*"], 110.95, 0.5],
      ["C5'", "C4'", "O4'", 109.47, 0.5],
      ["C5'", ["O5'","O5*"], "P", 120.73, 0.5],
      [["O5'","O5*"], "P", "OP1", 108.57, 0.5],
      [["O5'","O5*"], "P", "OP2", 107.42, 0.5],
      ["OP1", "P", "OP2", 119.86, 0.5],
      ["N9", "C1'", "O4'", 108.85, 0.5]
    ],
  }
};
function getAtomByName(resAtoms, nmArr){
  const set = new Set(nmArr.map(n=>n.toUpperCase()));
  for(const a of resAtoms){
    const nm = a.name.trim().toUpperCase();
    if(set.has(nm)) return a;
  }
  return null;
}
function _pick(resAtoms, nm){
  return Array.isArray(nm) ? getAtomByName(resAtoms, nm) : getAtomByName(resAtoms, [nm]);
}
function bondLen(resAtoms, aName, bName){
  const a = _pick(resAtoms, aName), b = _pick(resAtoms, bName);
  if(!a||!b) return null;
  return distance(a,b);
}
function angleDeg(resAtoms, aName, bName, cName){
  const a = _pick(resAtoms, aName);
  const b = _pick(resAtoms, bName);
  const c = _pick(resAtoms, cName);
  if(!a||!b||!c) return null;
  const ba = vnorm(vsub({x:a.x,y:a.y,z:a.z},{x:b.x,y:b.y,z:b.z}));
  const bc = vnorm(vsub({x:c.x,y:c.y,z:c.z},{x:b.x,y:b.y,z:b.z}));
  const cos = Math.min(1, Math.max(-1, vdot(ba,bc)));
  return (Math.acos(cos)*180/Math.PI);
}
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
  }else if(kind==="o5tors"){
    const ax = o5pAxisOf(resAtoms); if(!ax) return false;
    const subset = subsetByNames(resAtoms, PHOS_AROUND_O5P);
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
  const kinds = ["c5tors","o5tors","glyco"];
  let deltas = [2.0, 1.0, 0.5];
  for(const d of deltas){
    for(const rs of seqs){
      const resAtoms = atomsOf(atoms, chain, rs);
      let E0 = residueEnergy(resAtoms);
      let bestGain = 0, bestKind=null, bestSign=0;
      for(const kind of kinds){
        for(const sign of [+1,-1]){
          const snap = cloneCoords(resAtoms);
          const ok = trySmallRotation(resAtoms, kind, sign*d);
          if(!ok){ restoreCoords(snap); continue; }
          const E1 = residueEnergy(resAtoms);
          const gain = E0 - E1;
          if(gain > bestGain){ bestGain=gain; bestKind=kind; bestSign=sign; }
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
function spAtom(a){
  const nm = a.name.trim().toUpperCase();
  return SUGAR_NAMES.has(nm) || PHOS_NAMES.has(nm);
}
function vec(a){ return [a.x,a.y,a.z]; }
function add(u,v){ return [u[0]+v[0],u[1]+v[1],u[2]+v[2]]; }
function sub(u,v){ return [u[0]-v[0],u[1]-v[1],u[2]-v[2]]; }
function mul(u,s){ return [u[0]*s,u[1]*s,u[2]*s]; }
function dot(u,v){ return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]; }
function norm(u){ return Math.hypot(u[0],u[1],u[2]); }
const DIHEDRALS_IDEAL = {
  DNA: { alpha:-69.06402495, beta:-174.677887,  gamma:48.22070609, delta:141.1206416,  epsilon:-175.5731213, zeta:-95.64193751 },
  RNA: { alpha:-79.21083941, beta: 172.6145615, gamma:71.36607816, delta:110.001615,   epsilon: 179.7033593, zeta:-44.1887818  }
};
function periodicDeltaRad(target, current){
  let d = target - current;
  while (d >  Math.PI) d -= 2*Math.PI;
  while (d < -Math.PI) d += 2*Math.PI;
  return d;
}
function dihedralRad4(A,B,C,D){
  const b1 = [B.x-A.x, B.y-A.y, B.z-A.z];
  const b2 = [C.x-B.x, C.y-B.y, C.z-B.z];
  const b3 = [D.x-C.x, D.y-C.y, D.z-C.z];
  function unit(v){ const L=Math.hypot(v[0],v[1],v[2])||1; return [v[0]/L,v[1]/L,v[2]/L]; }
  function cross(u,v){ return [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]; }
  const n1 = unit(cross(b1,b2)), n2 = unit(cross(b2,b3));
  const m1 = cross(n1, unit(b2));
  const x = n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
  const y = m1[0]*n2[0]+m1[1]*n2[1]+m1[2]*n2[2];
  return Math.atan2(y, x);
}
function geometryBackboneEnergyGrad(atoms){
  const ideals = IDEALS_DATA[_CURRENT_NA_TYPE] || IDEALS_DATA.DNA;
  const N = atoms.length;
  const g = Array(N).fill(0).map(()=>[0,0,0]);
  let E=0;
  const idxMap = new Map(atoms.map((a,i)=>[a,i]));
  const grouped = _gbcrCache.get(atoms) || groupByChainRes(atoms);
  _gbcrCache.set(atoms, grouped);
  for (const key of grouped.keys()){
    const [chain,res] = key.split('|');
    const chainID = chain, resSeq = parseInt(res,10);
    const resAtoms = atomsOf(atoms, chainID, resSeq);
    for (const [A,B,d0,k] of ideals.BONDS){
      const a = _pick(resAtoms, A); const b = _pick(resAtoms, B);
      if(!a || !b) continue;
      const ia = idxMap.get(a), ib = idxMap.get(b);
      const ra = vec(a), rb = vec(b);
      const rij = sub(ra, rb);
      const d = norm(rij); if(d===0) continue;
      const diff = d - d0;
      E += (k||1.0)*diff*diff;
      const coeff = 2*(k||1.0)*diff/d;
      const gi = mul(rij, coeff), gj = mul(rij,-coeff);
      if (spAtom(a)) g[ia] = add(g[ia], gi);
      if (spAtom(b)) g[ib] = add(g[ib], gj);
    }
  }
  const eps = 1e-4;
  function angVal(a,b,c){
    const ba = sub(vec(a),vec(b)), bc = sub(vec(c),vec(b));
    const na = norm(ba), nc=norm(bc); if(na===0||nc===0) return Math.PI;
    let cosv = dot(ba,bc)/(na*nc); cosv = Math.max(-1,Math.min(1,cosv));
    return Math.acos(cosv); 
  }
  for (const key of grouped.keys()){
    const [chain,res] = key.split('|');
    const chainID = chain, resSeq = parseInt(res,10);
    const resAtoms = atomsOf(atoms, chainID, resSeq);
    for (const [A,B,C,a0deg,kw] of ideals.ANGLES){
      const Aat=_pick(resAtoms,A), Bat=_pick(resAtoms,B), Cat=_pick(resAtoms,C);
      if(!Aat||!Bat||!Cat) continue;
      const a0 = (a0deg*Math.PI/180);
      const theta = angVal(Aat,Bat,Cat);
      const diff = theta - a0; const w = (kw||0.5);
      E += w*diff*diff;

      const ids = [[Aat,'A'],[Bat,'B'],[Cat,'C']];
      for(const [atom,_label] of ids){
        if(!spAtom(atom)) continue;
        const ia = idxMap.get(atom);
        for(let c=0;c<3;c++){
          const keyc = c===0?'x':c===1?'y':'z';
          const bk = atom[keyc];
          atom[keyc] = bk + eps; const Ep = (function(){ const th=angVal(Aat,Bat,Cat); const df=th-a0; return w*df*df; })();
          atom[keyc] = bk - eps; const Em = (function(){ const th=angVal(Aat,Bat,Cat); const df=th-a0; return w*df*df; })();
          atom[keyc] = bk;
          g[ia][c] += (Ep-Em)/(2*eps);
        }
      }
    }
  }
  {
    const idealsT = DIHEDRALS_IDEAL[_CURRENT_NA_TYPE] || DIHEDRALS_IDEAL.DNA;
    const w_torsion = 0.12;
    const epsT = 1e-4;
    for (const key of grouped.keys()){
      const [chain,resStr] = key.split('|');
      const chainID = chain, resSeq = parseInt(resStr,10);
      const cur  = atomsOf(atoms, chainID, resSeq);
      const prev = atomsOf(atoms, chainID, resSeq-1);
      const next = atomsOf(atoms, chainID, resSeq+1);
      const pick = (arr, nm) => _pick(arr, nm);
      const sets = [
        ['alpha',   (prev&&cur) ? [ pick(prev,["O3'","O3*"]), pick(cur,["P"]), pick(cur,["O5'","O5*"]), pick(cur,["C5'"]) ] : null, idealsT.alpha ],
        ['beta',                [ pick(cur,["P"]), pick(cur,["O5'","O5*"]), pick(cur,["C5'"]), pick(cur,["C4'"]) ]                 , idealsT.beta  ],
        ['gamma',               [ pick(cur,["O5'","O5*"]), pick(cur,["C5'"]), pick(cur,["C4'"]), pick(cur,["C3'"]) ]                 , idealsT.gamma ],
        ['delta',               [ pick(cur,["C5'"]), pick(cur,["C4'"]), pick(cur,["C3'"]), pick(cur,["O3'","O3*"]) ]                 , idealsT.delta ],
        ['epsilon', (next) ?    [ pick(cur,["C4'"]), pick(cur,["C3'"]), pick(cur,["O3'","O3*"]), pick(next,["P"]) ]                 : null, idealsT.epsilon ],
        ['zeta',    (next) ?    [ pick(cur,["C3'"]), pick(cur,["O3'","O3*"]), pick(next,["P"]), pick(next,["O5'","O5*"]) ]          : null, idealsT.zeta   ],
      ];
      for (const ent of sets){
        if (!ent || !ent[1]) continue;
        const [name, arr4, tgtDeg] = ent;
        if (!arr4 || arr4.some(a=>!a)) continue;
        const [A,B,C,D] = arr4;
        const curRad = dihedralRad4(A,B,C,D);
        const tgtRad = tgtDeg * Math.PI/180.0;
        const dphi = periodicDeltaRad(tgtRad, curRad);
        E += w_torsion * dphi*dphi;
        const pts = [A,B,C,D];
        for (const P of pts){
          if(!spAtom(P)) continue;
          const i = idxMap.get(P);
          for (let c=0;c<3;c++){
            const keyc = c===0?'x':c===1?'y':'z';
            const bk = P[keyc];
            P[keyc] = bk + epsT;
            const Ep = (function(){
              const val=dihedralRad4(A,B,C,D); const df=periodicDeltaRad(tgtRad,val);
              return w_torsion*df*df;
            })();
            P[keyc] = bk - epsT;
            const Em = (function(){
              const val=dihedralRad4(A,B,C,D); const df=periodicDeltaRad(tgtRad,val);
              return w_torsion*df*df;
            })();
            P[keyc] = bk;
            g[i][c] += (Ep-Em)/(2*epsT);
          }
        }
      }
    }
  }
  return {E, grad:g};
}
function torsionalRepulsionFinal(atoms, {
  cutOP = 2.6, 
  cutPC = 3.0,
  maxIter = 12,
  stepDegTors = 2.0,
  stepDegAngle = 1.0 
} = {}){
  const grouped = _gbcrCache.get(atoms) || groupByChainRes(atoms);
  _gbcrCache.set(atoms, grouped);
  const chains = new Set(atoms.map(a=>a.chainID));
  for(const ch of chains){
    const seqs = residueSeqListForChain(atoms, ch);
    for(const rs of seqs){
      const res = atomsOf(atoms, ch, rs);
      const O5 = _pick(res, ["O5'","O5*"]);
      const P  = _pick(res, ["P"]);
      const OP1= _pick(res, ["OP1","O1P"]);
      const OP2= _pick(res, ["OP2","O2P"]);
      const C2p= _pick(res, ["C2'","C2*"]);
      const C3p= _pick(res, ["C3'","C3*"]);
      const O3p= _pick(res, ["O3'","O3*"]);
      const O2p= _pick(res, ["O2'","O2*"]); 
      if(!P || !O5 || !(OP1&&OP2)) continue;
      const o5pAxis = vnorm(vsub(P, O5));
      const originP = {x:P.x,y:P.y,z:P.z};
      function tryTorsionOnce(opAtom, targets, sign){
        const rad = (sign * stepDegTors) * Math.PI/180;
        const snap=[{a:opAtom,x:opAtom.x,y:opAtom.y,z:opAtom.z}];
        if (!rotateSubsetAroundAxis_radialOut([opAtom], originP, o5pAxis, rad)) return false;
        const ok = targets.every(t=>{
          if(!t) return true;
          return distance(opAtom, t) >= (t === C2p || t === C3p ? cutPC : cutOP);
        });
        if(!ok){
          for(const s of snap){ s.a.x=s.x; s.a.y=s.y; s.a.z=s.z; }
        }
        return ok;
      }
      function tryAngleOnce(opAtom, target, openToward='away'){
        if(!target) return true;
        const u = vnorm(vsub(opAtom, P));
        const v = vnorm(vsub(target, P));
        let axis = vcross(u, v);
        const nrm = Math.hypot(axis.x,axis.y,axis.z);
        if(nrm===0) return false;
        axis = vscale(axis, 1/nrm);
        const before = distance(opAtom, target);
        const rad1 = stepDegAngle * Math.PI/180;
        const snap=[{a:opAtom,x:opAtom.x,y:opAtom.y,z:opAtom.z}];
        if (!rotateSubsetAroundAxis_radialOut([opAtom], originP, axis, rad1)) {
          if (!rotateSubsetAroundAxis_radialOut([opAtom], originP, axis, -rad1)) {
            return false;
          }
        }
        const after1 = distance(opAtom, target);
        if(!(after1 > before)){
          for(const s of snap){ s.a.x=s.x; s.a.y=s.y; s.a.z=s.z; }
          rotateSubsetAroundAxis_radialOut([opAtom], originP, axis, -rad1);
          const after2 = distance(opAtom, target);
          if(!(after2 > before)){
            for(const s of snap){ s.a.x=s.x; s.a.y=s.y; s.a.z=s.z; }
            return false;
          }
        }
        return true;
      }
      for(let it=0; it<maxIter; it++){
        let improved = false;
        const neighborsOP = [C3p, O3p, C2p, O2p].filter(Boolean);
        const neighborsPC = [C2p, C3p].filter(Boolean);
        let needOP1 = neighborsOP.some(t => distance(OP1,t) < cutOP);
        if(needOP1){
          if(!neighborsOP.every(t => distance(OP1,t) >= cutOP)){
            improved = tryTorsionOnce(OP1, neighborsOP, +1) || tryTorsionOnce(OP1, neighborsOP, -1) || improved;
          }
          for(const t of neighborsOP){
            if(distance(OP1,t) < cutOP){
              improved = tryAngleOnce(OP1, t) || improved;
            }
          }
        }
        let needOP2 = neighborsOP.some(t => distance(OP2,t) < cutOP);
        if(needOP2){
          if(!neighborsOP.every(t => distance(OP2,t) >= cutOP)){
            improved = tryTorsionOnce(OP2, neighborsOP, +1) || tryTorsionOnce(OP2, neighborsOP, -1) || improved;
          }
          for(const t of neighborsOP){
            if(distance(OP2,t) < cutOP){
              improved = tryAngleOnce(OP2, t) || improved;
            }
          }
        }
        for(const t of neighborsPC){
          if(distance(P, t) < cutPC){
            improved = tryAngleOnce(OP1, t) || tryAngleOnce(OP2, t) || improved;
          }
        }
        if(!improved) break;
      }
    }
  }
}
function minimizeBackbone(
  atoms,
  { maxIter = 800, gtol = 1e-6, alpha0 = 0.20 } = {}
){
  let { E, grad } = geometryBackboneEnergyGrad(atoms);
  let step = alpha0;
  for (let it = 0; it < maxIter; it++) {
    let gnorm2 = 0;
    for (let i = 0; i < atoms.length; i++) {
      if (!spAtom(atoms[i])) continue;
      const gx = grad[i][0], gy = grad[i][1], gz = grad[i][2];
      gnorm2 += gx*gx + gy*gy + gz*gz;
    }
    const gnorm = Math.sqrt(gnorm2);
    if (gnorm < gtol) break;
    const dir = grad.map(v => [-v[0], -v[1], -v[2]]);
    const saved = atoms.map(a => [a.x, a.y, a.z]);
    const Eold  = E;
    let ok = false, alpha = step;
    for (let ls = 0; ls < 20; ls++) {
      for (let i = 0; i < atoms.length; i++) {
        if (!spAtom(atoms[i])) continue;
        atoms[i].x = saved[i][0] + alpha * dir[i][0];
        atoms[i].y = saved[i][1] + alpha * dir[i][1];
        atoms[i].z = saved[i][2] + alpha * dir[i][2];
      }
      const r = geometryBackboneEnergyGrad(atoms);
      if (r.E <= Eold - 1e-4 * alpha * gnorm2) { 
        E = r.E; grad = r.grad; ok = true; break;
      }
      alpha *= 0.5;
    }
    if (!ok) {
      for (let i = 0; i < atoms.length; i++) {
        atoms[i].x = saved[i][0];
        atoms[i].y = saved[i][1];
        atoms[i].z = saved[i][2];
      }
      step *= 0.5;
      if (step < 1e-6) break;
    }
  }
}
function finalContactSafeNudges(atoms, {
  cutOP = 2.6, 
  cutPC = 3.0, 
  passPerIter = 1,    
  stepDegTors = 0.8,  
  stepDegAngle = 0.4  
} = {}) {
  const grouped = _gbcrCache.get(atoms) || groupByChainRes(atoms);
  _gbcrCache.set(atoms, grouped);
  const chains = new Set(atoms.map(a=>a.chainID));
  for (let pass = 0; pass < passPerIter; pass++) {
    for (const ch of chains) {
      const seqs = residueSeqListForChain(atoms, ch);
      for (const rs of seqs) {
        const res = atomsOf(atoms, ch, rs);
        const O5 = _pick(res, ["O5'","O5*"]);
        const P  = _pick(res, ["P"]);
        const OP1= _pick(res, ["OP1","O1P"]);
        const OP2= _pick(res, ["OP2","O2P"]);
        const C2p= _pick(res, ["C2'","C2*"]);
        const C3p= _pick(res, ["C3'","C3*"]);
        const O3p= _pick(res, ["O3'","O3*"]);
        const O2p= _pick(res, ["O2'","O2*"]); 
        if (!P || !O5 || !(OP1 && OP2)) continue;
        const originP = {x:P.x, y:P.y, z:P.z};
        const u_o5p   = vnorm(vsub(P, O5));
        function nudgeTorsion(opAtom, sign){
          const rad = (sign * stepDegTors) * Math.PI/180;
          rotateSubsetAroundAxis_radialOut([opAtom], originP, u_o5p, rad);
        }
        function nudgeTiltAway(opAtom, target){
          if (!target) return;
          const u = vnorm(vsub(opAtom, P));
          const v = vnorm(vsub(target, P));
          let axis = vcross(u, v);
          const nrm = Math.hypot(axis.x, axis.y, axis.z);
          if (nrm === 0) return;
          axis = vscale(axis, 1/nrm);
          const before = distance(opAtom, target);
          const rad = stepDegAngle * Math.PI/180;
          const snap = {x:opAtom.x,y:opAtom.y,z:opAtom.z};
          if (!rotateSubsetAroundAxis_radialOut([opAtom], originP, axis, rad)) {
            rotateSubsetAroundAxis_radialOut([opAtom], originP, axis, -rad);
            return;
          }
          const after1 = distance(opAtom, target);
          if (!(after1 > before)) {
            opAtom.x = snap.x; opAtom.y = snap.y; opAtom.z = snap.z;
            rotateSubsetAroundAxis_radialOut([opAtom], originP, axis, -rad);
            const after2 = distance(opAtom, target);
            if (!(after2 > before)) {
              opAtom.x = snap.x; opAtom.y = snap.y; opAtom.z = snap.z;
            }
          }
        }
        const neighborsOP = [C3p, O3p, C2p, O2p].filter(Boolean);
        const neighborsPC = [C2p, C3p].filter(Boolean);
        for (const op of [OP1, OP2]) {
          let tooClose = neighborsOP.some(t => distance(op, t) < cutOP);
          if (tooClose) {
            const beforeMin = Math.min(...neighborsOP.map(t => distance(op, t)));
            nudgeTorsion(op, +1);
            let afterMin = Math.min(...neighborsOP.map(t => distance(op, t)));
            if (!(afterMin > beforeMin)) {
              nudgeTorsion(op, -1); nudgeTorsion(op, -1);
              afterMin = Math.min(...neighborsOP.map(t => distance(op, t)));
              if (!(afterMin > beforeMin)) {
                nudgeTorsion(op, +1);
              }
            }
          }
          for (const t of neighborsOP) {
            if (distance(op, t) < cutOP) nudgeTiltAway(op, t);
          }
        }
        for (const t of neighborsPC) {
          if (distance(P, t) < cutPC) {
            nudgeTiltAway(OP1, t);
            nudgeTiltAway(OP2, t);
          }
        }
      }
    }
  }
}
function enforceTetrahedralAroundP(atoms, {
  targetDeg = 109.5,
  stepDeg   = 0.5,
  iters     = 2,
  includeOP1 = true,
  includeOP2 = true
} = {}) {
  function angleAtP(P, X, Y){
    const u = vnorm(vsub(X, P));
    const v = vnorm(vsub(Y, P));
    let c = vdot(u, v);
    c = Math.max(-1, Math.min(1, c));
    return Math.acos(c) * 180/Math.PI;
  }
  function nudgeAngleToward(P, OP, other, stepDeg){
    if(!P || !OP || !other) return;
    const before = angleAtP(P, OP, other);
    const tgt = targetDeg;
    if (!isFinite(before)) return;
    const u = vnorm(vsub(OP, P));
    const v = vnorm(vsub(other, P));
    let axis = vcross(u, v);
    const nrm = Math.hypot(axis.x, axis.y, axis.z);
    if (nrm === 0) return;
    axis = vscale(axis, 1/nrm);
    const originP = {x:P.x, y:P.y, z:P.z};
    const snap = {x:OP.x, y:OP.y, z:OP.z};
    if (!rotateSubsetAroundAxis_radialOut([OP], originP, axis, ( stepDeg*Math.PI/180))) {
      OP.x = snap.x; OP.y = snap.y; OP.z = snap.z;
      if (!rotateSubsetAroundAxis_radialOut([OP], originP, axis, (-stepDeg*Math.PI/180))) {
        OP.x = snap.x; OP.y = snap.y; OP.z = snap.z;
        return;
      }
    }
    let afterPlus = angleAtP(P, OP, other);
    let plusGain = Math.abs(tgt - before) - Math.abs(tgt - afterPlus);
    OP.x = snap.x; OP.y = snap.y; OP.z = snap.z;
    rotateSubsetAroundAxis_radialOut([OP], originP, axis, (-stepDeg*Math.PI/180));
    let afterMinus = angleAtP(P, OP, other);
    let minusGain = Math.abs(tgt - before) - Math.abs(tgt - afterMinus);
    if (plusGain <= 0 && minusGain <= 0){
      OP.x = snap.x; OP.y = snap.y; OP.z = snap.z;
      return;
    }
    if (plusGain >= minusGain){
      OP.x = snap.x; OP.y = snap.y; OP.z = snap.z;
      rotateSubsetAroundAxis_radialOut([OP], originP, axis, ( stepDeg*Math.PI/180));
    }
  }
  const grouped = _gbcrCache.get(atoms) || groupByChainRes(atoms);
  _gbcrCache.set(atoms, grouped);
  const chains = new Set(atoms.map(a=>a.chainID));
  for (const ch of chains){
    const { pairs } = buildPairsForChain(atoms, ch);
    for (const [prevSeq, currSeq] of pairs){
      const prev = atomsOf(atoms, ch, prevSeq);
      const curr = atomsOf(atoms, ch, currSeq);
      const P   = _pick(curr, ["P"]);
      const O5  = _pick(curr, ["O5'","O5*"]);
      const OP1 = _pick(curr, ["OP1","O1P"]);
      const OP2 = _pick(curr, ["OP2","O2P"]);
      const O3p = _pick(prev, ["O3'","O3*"]);
      if (!P || !O5 || !(OP1 && OP2) || !O3p) continue;
      for (let it=0; it<iters; it++){
        if (includeOP1){
          nudgeAngleToward(P, OP1, O5,  stepDeg);
          nudgeAngleToward(P, OP1, OP2, stepDeg);
          nudgeAngleToward(P, OP1, O3p, stepDeg);
        }
        if (includeOP2){
          nudgeAngleToward(P, OP2, O5,  stepDeg);
          nudgeAngleToward(P, OP2, OP1, stepDeg);
          nudgeAngleToward(P, OP2, O3p, stepDeg);
        }
      }
    }
  }
}
function geometryBackboneFinalEnergyGrad(atoms, {wAngle=0.5, wPO=2.0}={}){
  const ideals = IDEALS_DATA[_CURRENT_NA_TYPE] || IDEALS_DATA.DNA;
  const N = atoms.length;
  const g = Array(N).fill(0).map(()=>[0,0,0]);
  let E=0;
  const idxMap = new Map(atoms.map((a,i)=>[a,i]));
  const grouped = _gbcrCache.get(atoms) || groupByChainRes(atoms);
  _gbcrCache.set(atoms, grouped);
  function addGradBond(a,b,d0,k=1.0){
    const ia = idxMap.get(a), ib = idxMap.get(b);
    const ra = vec(a), rb = vec(b);
    const rij = sub(ra, rb);
    const d = norm(rij); if(d===0) return;
    const diff = d - d0;
    E += k*diff*diff;
    const coeff = 2*k*diff/d;
    const gi = mul(rij, coeff), gj = mul(rij,-coeff);
    if (spAtom(a)) g[ia] = add(g[ia], gi);
    if (spAtom(b)) g[ib] = add(g[ib], gj);
  }
  function angRad(a,b,c){
    const ba = sub(vec(a),vec(b)), bc = sub(vec(c),vec(b));
    const na = norm(ba), nc=norm(bc); if(na===0||nc===0) return Math.PI;
    let cosv = dot(ba,bc)/(na*nc); cosv = Math.max(-1,Math.min(1,cosv));
    return Math.acos(cosv);
  }
  const eps = 1e-4;
  for (const key of grouped.keys()){
    const [chain,resS] = key.split('|');
    const chainID = chain, resSeq = parseInt(resS,10);
    const resAtoms = atomsOf(atoms, chainID, resSeq);
    for (const [A,B,d0,k] of ideals.BONDS){
      const a = _pick(resAtoms, A); const b = _pick(resAtoms, B);
      if(!a || !b) continue;
      addGradBond(a,b,d0,(k||1.0));
    }
  }
  for (const key of grouped.keys()){
    const [chain,resS] = key.split('|');
    const chainID = chain, resSeq = parseInt(resS,10);
    const resAtoms = atomsOf(atoms, chainID, resSeq);
    for (const [A,B,C,a0deg,_kw] of ideals.ANGLES){
      const Aat=_pick(resAtoms,A), Bat=_pick(resAtoms,B), Cat=_pick(resAtoms,C);
      if(!Aat||!Bat||!Cat) continue;
      const a0 = (a0deg*Math.PI/180);
      const theta = angRad(Aat,Bat,Cat);
      const diff = theta - a0;
      E += wAngle*diff*diff;
      const ids = [Aat,Bat,Cat];
      for(const atom of ids){
        if(!spAtom(atom)) continue;
        const ia = idxMap.get(atom);
        for(let c=0;c<3;c++){
          const keyc = c===0?'x':c===1?'y':'z';
          const bk = atom[keyc];
          atom[keyc] = bk + eps; const Ep = (function(){ const th=angRad(Aat,Bat,Cat); const df=th-a0; return wAngle*df*df; })();
          atom[keyc] = bk - eps; const Em = (function(){ const th=angRad(Aat,Bat,Cat); const df=th-a0; return wAngle*df*df; })();
          atom[keyc] = bk;
          g[ia][c] += (Ep-Em)/(2*eps);
        }
      }
    }
  }
  const chains = new Set(atoms.map(a=>a.chainID));
  for(const ch of chains){
    const {pairs} = buildPairsForChain(atoms, ch);
    for(const [prevSeq,currSeq] of pairs){
      const prev = atomsOf(atoms,ch,prevSeq);
      const curr = atomsOf(atoms,ch,currSeq);
      const O3p = _pick(prev, ["O3'","O3*"]);
      const P   = _pick(curr, ["P"]);
      if(!O3p||!P) continue;
      addGradBond(O3p, P, 1.60, wPO);
    }
  }
  return {E, grad:g};
}
function minimizeBackboneFinalLSQ(atoms, {maxIter=800, gtol=1e-6, alpha0=0.18, wAngle=0.5, wPO=2.0}={}) {
  let {E, grad} = geometryBackboneFinalEnergyGrad(atoms, {wAngle, wPO});
  let step = alpha0;
  for(let it=0; it<maxIter; it++){
    finalContactSafeNudges(atoms, {
      cutOP: 3.2,
      cutPC: 3.5,
      passPerIter: 1,
      stepDegTors: 0.8,
      stepDegAngle: 0.4
    });
    enforceTetrahedralAroundP(atoms, {
      targetDeg: 109.5,
      stepDeg:   0.4, 
      iters:     1,  
      includeOP1: true,
      includeOP2: true
    });
    let gnorm=0;
    for(let i=0;i<atoms.length;i++){
      if(!spAtom(atoms[i])) continue;
      gnorm += dot(grad[i],grad[i]);
    }
    gnorm = Math.sqrt(gnorm);
    if (gnorm < gtol) break;
    const dir = grad.map(v=>[-v[0],-v[1],-v[2]]);
    const saved = atoms.map(a=>[a.x,a.y,a.z]);
    const Eold = E;
    let ok=false, alpha=step;
    for(let ls=0; ls<20; ls++){
      for(let i=0;i<atoms.length;i++){
        if(!spAtom(atoms[i])) continue;
        atoms[i].x = saved[i][0] + alpha*dir[i][0];
        atoms[i].y = saved[i][1] + alpha*dir[i][1];
        atoms[i].z = saved[i][2] + alpha*dir[i][2];
      }
      const r = geometryBackboneFinalEnergyGrad(atoms, {wAngle, wPO});
      if (r.E <= Eold - 1e-4*alpha*gnorm*gnorm){ E=r.E; grad=r.grad; ok=true; break; }
      alpha *= 0.5;
    }
    if(!ok){
      for(let i=0;i<atoms.length;i++){ atoms[i].x=saved[i][0]; atoms[i].y=saved[i][1]; atoms[i].z=saved[i][2]; }
      step *= 0.5;
      if(step < 1e-6) break;
    }
    finalContactSafeNudges(atoms);
  }
}
function vdwRadiusOf(atom){
  const el = (atom.element && atom.element.trim()) ? atom.element.trim().toUpperCase()
            : atom.name.trim().slice(0,1).toUpperCase();
  switch(el){
    case 'H': return 1.20;
    case 'C': return 1.70;
    case 'N': return 1.55;
    case 'O': return 1.52;
    case 'P': return 1.80;
    case 'S': return 1.80;
    default:  return 1.70;
  }
}
function buildAdjacencyForResidue(resAtoms, ideals){
  const adj = new Map();
  function add(n1,n2){
    const a=(Array.isArray(n1)?n1:[n1]), b=(Array.isArray(n2)?n2:[n2]);
    for(const x of a){ for(const y of b){
      const X = x.toUpperCase(), Y = y.toUpperCase();
      if(!adj.has(X)) adj.set(X,new Set());
      if(!adj.has(Y)) adj.set(Y,new Set());
      adj.get(X).add(Y); adj.get(Y).add(X);
    }}
  }
  for(const [A,B] of ideals.BONDS){ add(A,B); }
  return adj;
}
function makePO3LinkSkipSet(atoms){
  const skip = new Set();
  const chains = new Set(atoms.map(a=>a.chainID));
  for(const ch of chains){
    const { pairs } = buildPairsForChain(atoms, ch);
    for (const [prevSeq, currSeq] of pairs){
      const key = `${ch}|${prevSeq}|O3-${ch}|${currSeq}|P`;
      skip.add(key);
    }
  }
  return skip;
}
function isBondedWithinResidue(a, b, ideals, grouped){
  if (a.chainID!==b.chainID || a.resSeq!==b.resSeq) return false;
  const key = `${a.chainID}|${a.resSeq}`;
  const resAtoms = grouped.get(key) || [];
  const adj = buildAdjacencyForResidue(resAtoms, ideals);
  const A = a.name.trim().toUpperCase(), B = b.name.trim().toUpperCase();
  return (adj.get(A)?.has(B) || false);
}
function resolveClashesLocal(atoms, {
  cutoff = 4.5,
  scaleS = 0.90,
  betaNeighbor = 0.25,
  iters = 10,
  maxShiftPerIter = 0.25
} = {}){
  const ideals = IDEALS_DATA[_CURRENT_NA_TYPE] || IDEALS_DATA.DNA;
  const idxMap = new Map(atoms.map((a,i)=>[a,i]));
  const grouped = _gbcrCache.get(atoms) || groupByChainRes(atoms);
  _gbcrCache.set(atoms, grouped);
  const skipPO = makePO3LinkSkipSet(atoms);
  function resKeyOf(a){ return `${a.chainID}|${a.resSeq}`; }
  function adjFor(resKey){
    const resAtoms = grouped.get(resKey) || [];
    return buildAdjacencyForResidue(resAtoms, ideals);
  }
  for(let it=0; it<iters; it++){
    const disp = Array(atoms.length).fill(0).map(()=>({x:0,y:0,z:0}));
    for(let i=0; i<atoms.length; i++){
      const ai = atoms[i]; if(!spAtom(ai)) continue;
      for(let j=i+1; j<atoms.length; j++){
        const aj = atoms[j]; if(!spAtom(aj)) continue;
        const dx = aj.x-ai.x, dy=aj.y-ai.y, dz=aj.z-ai.z;
        const r  = Math.hypot(dx,dy,dz);
        if (r===0 || r>cutoff) continue;
        if (isBondedWithinResidue(ai, aj, ideals, grouped)) continue;
        const nmI = ai.name.trim().toUpperCase(), nmJ = aj.name.trim().toUpperCase();
        const k1 = `${ai.chainID}|${ai.resSeq}|${nmI.startsWith('O3')?'O3':nmI}-${aj.chainID}|${aj.resSeq}|${nmJ==='P'?'P':nmJ}`;
        const k2 = `${ai.chainID}|${ai.resSeq}|${nmI==='P'?'P':nmI}-${aj.chainID}|${aj.resSeq}|${nmJ.startsWith('O3')?'O3':nmJ}`;
        if (skipPO.has(k1) || skipPO.has(k2)) continue;
        const ri = vdwRadiusOf(ai), rj = vdwRadiusOf(aj);
        const rmin = scaleS * (ri + rj);
        if (r >= rmin) continue;
        let overlap = rmin - r;
        if (overlap > 2*maxShiftPerIter) overlap = 2*maxShiftPerIter;
        const ux = dx / r, uy = dy / r, uz = dz / r;
        const s  = 0.5 * overlap;
        const di = {x: -s*ux, y: -s*uy, z: -s*uz};
        const dj = {x: +s*ux, y: +s*uy, z: +s*uz};
        const ii = idxMap.get(ai), jj = idxMap.get(aj);
        disp[ii].x += di.x; disp[ii].y += di.y; disp[ii].z += di.z;
        disp[jj].x += dj.x; disp[jj].y += dj.y; disp[jj].z += dj.z;
      }
    }
    for(let k=0;k<atoms.length;k++){
      if(!spAtom(atoms[k])) continue;
      atoms[k].x += disp[k].x;
      atoms[k].y += disp[k].y;
      atoms[k].z += disp[k].z;
    }
  }
}
function runPOFitTwoPassesForChain(atoms, chain){
  const seqs = residueSeqListForChain(atoms, chain);
  if (seqs.length < 2) return;
  const { pairs } = buildPairsForChain(atoms, chain);
  for (let pass=0; pass<2; pass++){
    for (const [prevSeq, currSeq] of pairs){
      const prev = atomsOf(atoms, chain, prevSeq);
      const curr = atomsOf(atoms, chain, currSeq);
      const O3p = _pick(prev, ["O3'","O3*"]);
      const P   = _pick(curr, ["P"]);
      if (O3p && P) {
        const d = distance(O3p, P);
        TARGET_PO = (d && d > 2.6) ? Math.min(d - 0.4, 2.2) : 1.60;
      } else {
        TARGET_PO = 1.60;
      }
      const plan_g = bestAngle_Glyco(prev, curr);
      apply_Glyco(curr, plan_g);
      const plan_c5 = bestAngle_C5torsion(prev, curr);
      apply_C5torsion(curr, plan_c5);
      const plan_o5 = bestAngle_O5torsion(prev, curr);
      apply_O5torsion(curr, plan_o5);
    }
  }
}
function runPOFitAllChains(atoms){
  const chains = new Set(atoms.map(a=>a.chainID));
  for(const ch of chains) runPOFitTwoPassesForChain(atoms, ch);
}
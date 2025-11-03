const TARGET_PO = 1.60;
const SUGAR_NAMES = new Set([
  "C1'", "C2'", "C3'", "C4'", "O4'", "C5'", "O5'", "O3'", "O2'",
  "C1*", "C2*", "C3*", "C4*", "O4*", "C5*", "O5*", "O3*", "O2*"
]);
const PHOS_NAMES = new Set(["P","OP1","OP2","OP3","O1P","O2P","O3P"]);
const PHOS_GROUP = new Set(["P","OP1","OP2","OP3","O1P","O2P","O3P","O5'","O5*"]);
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
function findAtom(atoms, names){
  const set = new Set(names.map(n=>n.toUpperCase()));
  for(const a of atoms){
    const nm = a.name.trim().toUpperCase();
    if(set.has(nm)) return a;
  }
  return null;
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
  const axis   = vnorm(vsub({x:O5p.x,y:O5p.y,z:O5p.z}, origin)); // C5′→O5′
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
  let bestDeg = 0, bestCost = Infinity
  for(let deg=0; deg<360; deg+=3){
    const rad = deg*Math.PI/180;
    const Pp  = rotateAroundAxis({x:P.x,y:P.y,z:P.z}, ax.origin, ax.axis, rad);
    const d = distance(Pp, O3p);
    const cost = (d - TARGET_PO) * (d - TARGET_PO);
    if(cost < bestCost){ bestCost = cost; bestDeg = deg; }
  }
  let fine = bestDeg, fineCost = bestCost;
  for(let deg=bestDeg-3; deg<=bestDeg+3; deg+=0.2){
    let nd = (deg<0?deg+360:deg%360);
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




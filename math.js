function centroid(atoms) {
  let sx=0, sy=0, sz=0, n=atoms.length;
  for (const a of atoms) { sx+=a.x; sy+=a.y; sz+=a.z; }
  return (n? {x:sx/n, y:sy/n, z:sz/n} : {x:0,y:0,z:0});
}
function vsub(a,b){ return {x:a.x-b.x,y:a.y-b.y,z:a.z-b.z}; }
function vadd(a,b){ return {x:a.x+b.x,y:a.y+b.y,z:a.z+b.z}; }
function vdot(a,b){ return a.x*b.x+a.y*b.y+a.z*b.z; }
function vcross(a,b){ return {x:a.y*b.z-a.z*b.y, y:a.z*b.x-a.x*b.z, z:a.x*b.y-a.y*b.x}; }
function vscale(a,s){ return {x:a.x*s,y:a.y*s,z:a.z*s}; }
function vnorm(a){ const n=Math.hypot(a.x,a.y,a.z); return n? vscale(a,1/n) : {x:0,y:0,z:0}; }
function rotateAroundAxis(point, origin, axisUnit, rad){
  const po = vsub(point, origin);
  const c = Math.cos(rad), s = Math.sin(rad);
  const term1 = vscale(po, c);
  const term2 = vscale(vcross(axisUnit, po), s);
  const term3 = vscale(axisUnit, vdot(axisUnit, po)*(1-c));
  return vadd(origin, vadd(term1, vadd(term2, term3)));
}
function rotateAtomsAroundAxis(atoms, origin, axisUnit, rad){
  for(const a of atoms){
    const p = {x:a.x,y:a.y,z:a.z};
    const pr = rotateAroundAxis(p, origin, axisUnit, rad);
    a.x = pr.x; a.y = pr.y; a.z = pr.z;
  }
}
function rotateSubsetAroundAxis(atomsSubset, origin, axisUnit, rad){
  for(const a of atomsSubset){
    const p = {x:a.x,y:a.y,z:a.z};
    const pr = rotateAroundAxis(p, origin, axisUnit, rad);
    a.x = pr.x; a.y = pr.y; a.z = pr.z;
  }
}
function radialR(atom){ return Math.hypot(atom.x, atom.y); }
function rotateSubsetAroundAxis_radialOut(atomsSubset, origin, axisUnit, rad){
  const snap = atomsSubset.map(a=>({a,x:a.x,y:a.y,z:a.z,r0:radialR(a)}));
  rotateSubsetAroundAxis(atomsSubset, origin, axisUnit, rad);
  let ok = true;
  for (const s of snap){
    if (radialR(s.a) <= s.r0 - 1e-6) { ok=false; break; }
  }
  if (!ok){
    rotateSubsetAroundAxis(atomsSubset, origin, axisUnit, -2*rad);
    ok = true;
    for (const s of snap){
      if (radialR(s.a) <= s.r0 - 1e-6) { ok=false; break; }
    }
    if (!ok){ for (const s of snap){ s.a.x=s.x; s.a.y=s.y; s.a.z=s.z; } }
  }
  return ok;
}
function isOutwardOP(OP, P){
  const dx = OP.x - P.x, dy = OP.y - P.y;
  const rx = P.x,        ry = P.y;
  return (dx*rx + dy*ry) >= 0;
}
function rotateOP_withOutward(OP, origin, axisUnit, rad, P){
  const snap = {x:OP.x,y:OP.y,z:OP.z};
  rotateSubsetAroundAxis([OP], origin, axisUnit, rad);
  if(!isOutwardOP(OP, P)){
    OP.x=snap.x; OP.y=snap.y; OP.z=snap.z;
    rotateSubsetAroundAxis([OP], origin, axisUnit, -rad);
    if(!isOutwardOP(OP, P)){
      OP.x=snap.x; OP.y=snap.y; OP.z=snap.z;
      return false;
    }
  }
  return true;
}
function findAllAtoms(atoms, names){
  const want = new Set(names.map(s=>s.toUpperCase()));
  return atoms.filter(a => want.has(a.name.trim().toUpperCase()));
}
function meanC1pXYZ(layerAtoms){
  const c1s = findAllAtoms(layerAtoms, ["C1'", "C1*"]);
  if (!c1s.length) return null;
  let sx=0, sy=0, sz=0;
  for (const a of c1s){ sx+=a.x; sy+=a.y; sz+=a.z; }
  const n=c1s.length;
  return {x:sx/n, y:sy/n, z:sz/n};
}
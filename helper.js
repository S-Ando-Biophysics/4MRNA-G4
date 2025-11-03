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
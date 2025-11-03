function parseAtomLine(line) {
  const rec = line.slice(0,6).trim();
  if (!(rec === "ATOM" || rec === "HETATM")) return null;
  return {
    rec,
    serial: parseInt(line.slice(6,11)),
    name: line.slice(12,16),
    altLoc: line.slice(16,17),
    resName: line.slice(17,20),
    chainID: line.slice(21,22) || "A",
    resSeq: parseInt(line.slice(22,26)),
    iCode: line.slice(26,27),
    x: parseFloat(line.slice(30,38)),
    y: parseFloat(line.slice(38,46)),
    z: parseFloat(line.slice(46,54)),
    occ: line.slice(54,60).trim(),
    temp: line.slice(60,66).trim(),
    element: line.length >= 78 ? line.slice(76,78).trim() : "",
    charge: line.length >= 80 ? line.slice(78,80).trim() : "",
  };
}
function formatAtom(a) {
  const rec     = (a.rec || "ATOM").padEnd(6, " ");
  const serial  = a.serial.toString().padStart(5, " ");
  const name    = a.name.padStart(4, " ");
  const altLoc  = (a.altLoc || " ").slice(0,1);
  const resName = (a.resName || "   ").padStart(3, " ");
  const chainID = (a.chainID || "A").slice(0,1);
  const resSeq  = (a.resSeq || 1).toString().padStart(4, " ");
  const iCode   = (a.iCode || " ").slice(0,1);
  const x       = a.x.toFixed(3).toString().padStart(8, " ");
  const y       = a.y.toFixed(3).toString().padStart(8, " ");
  const z       = a.z.toFixed(3).toString().padStart(8, " ");
  const occ     = (a.occ && a.occ.length ? a.occ : "1.00").toString().padStart(6, " ");
  const temp    = (a.temp && a.temp.length ? a.temp : "0.00").toString().padStart(6, " ");
  const blanks67_76 = " ".repeat(10);
  const element = (a.element || "").toString().trim().padStart(2, " ");
  const charge  = (a.charge || "").toString().padStart(2, " ");
  return rec + serial + " " + name + altLoc + resName + " " + chainID + resSeq + iCode + "   " +
         x + y + z + occ + temp + blanks67_76 + element + charge;
}
function parsePDB(text) {
  const atoms = [];
  const carry = [];
  const lines = text.split(/\r?\n/);
  for (const line of lines) {
    const rec = line.slice(0,6).trim();
    if (rec === "ATOM" || rec === "HETATM") {
      const a = parseAtomLine(line); if (a) atoms.push(a);
    } else if (line.trim().length) {
      carry.push(line);
    }
  }
  return { atoms, carry };
}
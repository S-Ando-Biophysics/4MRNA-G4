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
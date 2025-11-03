(function(global){
  const crcTable = new Uint32Array(256);
  for(let n=0;n<256;n++){
    let c = n; for(let k=0;k<8;k++) c = (c & 1) ? (0xEDB88320 ^ (c >>> 1)) : (c >>> 1);
    crcTable[n] = c >>> 0;
  }
  function crc32(buf){ let c = 0 ^ (-1); for(let i=0;i<buf.length;i++) c = (c >>> 8) ^ crcTable[(c ^ buf[i]) & 0xFF]; return (c ^ (-1)) >>> 0; }
  const enc = new TextEncoder();
  const toU8 = (t) => enc.encode(String(t));
  const u16 = (v) => { const b = new Uint8Array(2); new DataView(b.buffer).setUint16(0,v,true); return b; };
  const u32 = (v) => { const b = new Uint8Array(4); new DataView(b.buffer).setUint32(0,v,true); return b; };
  function dosDateTime(){ return {time:0, date: (1<<5)|1}; }
  function concat(arrs){ let L=0; for(const a of arrs) L+=a.length; const out=new Uint8Array(L); let o=0; for(const a of arrs){ out.set(a,o); o+=a.length; } return out; }
  function lfHeader(nameBytes, dataBytes, crc){
    const {time,date}=dosDateTime();
    return concat([ u32(0x04034b50), u16(20), u16(0), u16(0), u16(time), u16(date),
                    u32(crc), u32(dataBytes.length), u32(dataBytes.length),
                    u16(nameBytes.length), u16(0) ]);
  }
  function cfHeader(nameBytes, dataBytes, crc, offset){
    const {time,date}=dosDateTime();
    return concat([ u32(0x02014b50), u16(20), u16(20), u16(0), u16(0), u16(time), u16(date),
                    u32(crc), u32(dataBytes.length), u32(dataBytes.length),
                    u16(nameBytes.length), u16(0), u16(0), u16(0), u16(0), u32(0), u32(offset) ]);
  }
  function eocd(count, size, offset){
    return concat([ u32(0x06054b50), u16(0), u16(0), u16(count), u16(count), u32(size), u32(offset), u16(0) ]);
  }
  function buildZip(files){
    const chunks=[], centers=[]; let offset=0;
    for(const f of files){
      const nameBytes=toU8(f.name);
      const dataBytes=f.text instanceof Uint8Array?f.text:toU8(f.text);
      const crc=crc32(dataBytes);
      const lh=lfHeader(nameBytes,dataBytes,crc);
      chunks.push(lh,nameBytes,dataBytes);
      const ch=cfHeader(nameBytes,dataBytes,crc,offset);
      centers.push(ch,nameBytes);
      offset += lh.length + nameBytes.length + dataBytes.length;
    }
    const central = concat(centers);
    const end = eocd(files.length, central.length, offset);
    const zipBytes = concat([...chunks, central, end]);
    return new Blob([zipBytes], {type:'application/zip'});
  }
  async function saveZip(files, zipname){
    const blob=buildZip(files);
    const url=URL.createObjectURL(blob);
    const a=document.createElement('a'); a.href=url; a.download=zipname||'bundle.zip';
    document.body.appendChild(a); a.click(); a.remove();
    URL.revokeObjectURL(url);
  }
  global.zipAndSave = async (files, name) => saveZip(files, name);
})(window);

(function(){
  const handBatchStatus = document.getElementById('handBatchStatus') || document.getElementById('status');
  const numLayersEl = document.getElementById('numLayers');
  const globalHandEl = document.getElementById('globalHand');
  const topologyEl = document.getElementById('topology');
  const naTypeEl = document.getElementById('naType');
  function getLayerPattern(N, hand, topology){
    const cw='cw', acw='acw';
    const isRight = (hand === 'right');
    const L = [];
    switch (topology) {
      case 'parallel':
        for (let i=0;i<N;i++) L.push(isRight ? cw : acw);
        break;
      case 'antiparallel-basket':
      case 'antiparallel-chair': {
        const start = isRight ? acw : cw;
        for (let i=0;i<N;i++) L.push((i%2===0)? start : (start===cw?acw:cw));
        break;
      }
      case 'hybrid':
        if (isRight) { for (let i=0;i<N;i++) L.push(cw); }
        else { for (let i=0;i<N;i++) L.push((i%2===0)? cw : acw); }
        break;
      default:
        for (let i=0;i<N;i++) L.push(cw);
    }
    return L;
  }
  function* tripleCombos(gaps){
    const total = Math.pow(3,gaps);
    for(let t=0;t<total;t++){
      const arr=new Array(gaps).fill(0);
      let x=t; for(let i=0;i<gaps;i++){ arr[i]=x%3; x=Math.floor(x/3); }
      yield arr;
    }
  }
  const Ltag = L => L.map(x => x==='acw' ? 'A' : 'C').join('');
  function Ttag(L,H,tripleIdxArr){
    if (!H.length) return '';
    const arr=[];
    for(let i=0;i<H.length;i++){
      const key = `${L[i]}->${L[i+1]}`;
      const hand = H[i];
      const triRaw = TWISTS_RAW[hand][key];
      const idx = (Array.isArray(tripleIdxArr) ? tripleIdxArr[i] : 1) || 0;
      arr.push(triRaw[idx].toFixed(2));
    }
    return arr.join('-');
  }
  function setSelectAndDispatch(sel, value){
    if (!sel) return;
    sel.value = value;
    sel.dispatchEvent(new Event('change'));
  }
  async function runAll(kind){
    try{
      const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value||'1')));
      const gaps = Math.max(0, N-1);
      const sels = document.querySelectorAll('.layerTemplate');
      const originalHand = globalHandEl.value;
      const originalTopo = topologyEl?.value;
      const topoList = ['parallel','antiparallel-basket','antiparallel-chair','hybrid'];
      const handVal  = (kind === 'left') ? 'left' : 'right';
      const H = new Array(gaps).fill(handVal);
      const allFiles = [];
      let gidx = 0; 
      setSelectAndDispatch(globalHandEl, handVal);
      for (const topology of topoList) {
        setSelectAndDispatch(topologyEl, topology);
        const L = getLayerPattern(N, handVal, topology);
        L.forEach((v,i)=> setSelectAndDispatch(sels[i], v));
        applyRiseAndTwistFromRules(L, H, null);
        for (const T of tripleCombos(gaps)){
          applyRiseAndTwistFromRules(L, H, T);
          try {
            const pdb = await generatePDB();
            const fname = `${String(++gidx).padStart(4,'0')}_${N}layers_${topology}_${Ltag(L)}_${Ttag(L,H,T)}.pdb`;
            allFiles.push({ name: fname, text: pdb });
            if (gidx % 50 === 0 && handBatchStatus) handBatchStatus.textContent = `Generating (${kind}): ${gidx}`;
            await new Promise(r=>setTimeout(r,0));
          } catch(e) {
            continue;
          }
        }
      }
      const zipName = (kind === 'right') ? '4MRNA-G4-Models-Right.zip' : '4MRNA-G4-Models-Left.zip';
      await window.zipAndSave(allFiles, zipName);
      if (handBatchStatus) handBatchStatus.textContent = `Done: saved ${allFiles.length} structures to ZIP.`;
      setSelectAndDispatch(globalHandEl, originalHand);
      if (topologyEl && originalTopo!=null) setSelectAndDispatch(topologyEl, originalTopo);
    }catch(e){
      const msg = 'Error: ' + (e?.message||e);
      if (handBatchStatus) handBatchStatus.textContent = msg;
      alert(msg);
    }
  }
  document.getElementById('dlAllRight').addEventListener('click', ()=>runAll('right'));
  document.getElementById('dlAllLeft').addEventListener('click',  ()=>runAll('left'));
})();
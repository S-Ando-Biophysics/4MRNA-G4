let viewer = null;
function ensureViewer() {
  if (!viewer) viewer = $3Dmol.createViewer(document.getElementById("viewer"), { backgroundColor: 'white' });
}
function showSpinner(){ const sp=document.getElementById('viewerSpinner'); if(sp) sp.style.display='flex'; }
function hideSpinner(){ const sp=document.getElementById('viewerSpinner'); if(sp) sp.style.display='none'; }
async function waitNextFrame(){ return new Promise(res=>requestAnimationFrame(()=>res())); }
function resetViewerNow(){
  ensureViewer();
  try { viewer.clear(); viewer.render(); } catch(_){}
}
function draw3D(pdbText, N) {
  ensureViewer();
  viewer.clear();
  try {
    viewer.addModel(pdbText, 'pdb');
  } catch (e) {
    const statusEl = document.getElementById('status');
    if (statusEl) statusEl.textContent = 'Error';
    hideSpinner();
    return;
  }
  const chainColors = { 'A': 0xff00ff, 'B': 0x00ffff, 'C': 0x00ff00, 'D': 0xffa500 };
  ['A','B','C','D'].forEach(ch => {
    const baseColor = chainColors[ch] || 0x808080;
    viewer.setStyle({ chain: ch },           { stick: { radius: 0.15, color: baseColor } });
    viewer.setStyle({ chain: ch, elem: 'N' },{ stick: { radius: 0.15, color: 0x0000ff } });
    viewer.setStyle({ chain: ch, elem: 'O' },{ stick: { radius: 0.15, color: 0xff0000 } });
  });
  viewer.zoomTo();
  viewer.render();
  setTimeout(() => {
    try { viewer.resize(); viewer.zoomTo(); viewer.render(); }
    catch(_){} finally { hideSpinner(); }
  }, 0);
}
window.addEventListener("resize", () => { if (viewer) { try{ viewer.resize(); viewer.render(); }catch(_){} } });
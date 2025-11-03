document.getElementById("downloadBtn").addEventListener("click", () => {
  if (!outText) return;
  const na = (naTypeEl?.value || 'DNA').toUpperCase();
  const topology = topologyEl?.value || 'parallel';
  const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value || "1")));
  const hand = globalHandEl.value || 'right';
  const fname = `G4_${na}_${N}layers_${hand}_${topology}.pdb`;
  const blob = new Blob([outText], { type: "chemical/x-pdb" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = fname;
  document.body.appendChild(a); a.click(); a.remove();
  URL.revokeObjectURL(url);
});
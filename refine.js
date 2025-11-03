(function () {
  const refineBtn = document.getElementById("refineBtn");
  if (!refineBtn) return;
  refineBtn.addEventListener("click", async () => {
    btnDl.disabled = true;
    try {
      outText = await generatePDB(); 
      const N = Math.max(1, Math.min(4, parseInt(numLayersEl.value || "1")));
      draw3D(outText, N);
      btnDl.disabled = false;
    } catch (_) {
    }
  });
})();


document.addEventListener('DOMContentLoaded', () => {
  const mrBtn = document.getElementById('mrToggleBtn');
  const batchSection = document.getElementById('batchDownloadSection');
  const progressWrap = document.getElementById('batchProgress');
  const bar = document.getElementById('batchProgressBar');
  const barText = document.getElementById('batchProgressText');
  const statusEl = document.getElementById('handBatchStatus') || document.getElementById('status');
  function resetProgress() {
    if (statusEl) statusEl.textContent = '';
    if (progressWrap) progressWrap.style.display = 'none';
    if (bar) {
      bar.style.width = '0%';
      bar.setAttribute('aria-valuenow', '0');
    }
    if (barText) barText.textContent = '0%';
  }
  mrBtn?.addEventListener('click', () => {
    const hidden = batchSection.style.display === 'none' || batchSection.style.display === '';
    batchSection.style.display = hidden ? 'block' : 'none';
    if (hidden) resetProgress();
  });
});
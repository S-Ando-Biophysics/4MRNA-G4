document.addEventListener('DOMContentLoaded', () => {
  const dlPhaserBtn = document.getElementById('dlPhaser');
  const PHASER_RAW_URLS = [
    'https://raw.githubusercontent.com/S-Ando-Biophysics/4MRNA-G4/main/Templates/4MRNA-G4.sh',
    'https://raw.githubusercontent.com/S-Ando-Biophysics/4MRNA-G4/refs/heads/main/Templates/4MRNA-G4.sh'
  ];
  async function downloadPhaserScript() {
    dlPhaserBtn.disabled = true;
    let text = null;
    for (const url of PHASER_RAW_URLS) {
      try {
        const resp = await fetch(url, { cache: 'no-store' });
        if (resp.ok) {
          text = await resp.text();
          break;
        }
      } catch (e) {
      }
    }
    if (text) {
      const blob = new Blob([text], { type: 'application/x-sh' });
      const link = document.createElement('a');
      link.href = URL.createObjectURL(blob);
      link.download = '4MRNA-G4.sh';
      document.body.appendChild(link);
      link.click();
      URL.revokeObjectURL(link.href);
      link.remove();
    } else {
      alert('Failed to download 4MRNA-G4.sh');
    }
    dlPhaserBtn.disabled = false;
  }
  dlPhaserBtn?.addEventListener('click', downloadPhaserScript);
});
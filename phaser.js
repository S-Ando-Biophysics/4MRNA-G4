document.addEventListener('DOMContentLoaded', () => {
  const dlPhaserOneBtn = document.getElementById('dlPhaserOne');
  const dlPhaserTwoBtn = document.getElementById('dlPhaserTwo');
  const PHASER_URLS_ONE = [
    'https://raw.githubusercontent.com/S-Ando-Biophysics/4MRNA-G4/main/Templates/4MRNA-G4-OneBlock.sh',
    'https://raw.githubusercontent.com/S-Ando-Biophysics/4MRNA-G4/refs/heads/main/Templates/4MRNA-G4-OneBlock.sh'
  ];
  const PHASER_URLS_TWO = [
    'https://raw.githubusercontent.com/S-Ando-Biophysics/4MRNA-G4/main/Templates/4MRNA-G4-TwoBlock.sh',
    'https://raw.githubusercontent.com/S-Ando-Biophysics/4MRNA-G4/refs/heads/main/Templates/4MRNA-G4-TwoBlock.sh'
  ];
  async function downloadPhaserScript(urlList, filename, btn) {
    btn.disabled = true;
    let text = null;
    for (const url of urlList) {
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
      link.download = filename;
      document.body.appendChild(link);
      link.click();
      URL.revokeObjectURL(link.href);
      link.remove();
    } else {
      alert(`Failed to download ${filename}`);
    }
    btn.disabled = false;
  }
  dlPhaserOneBtn?.addEventListener('click', () =>
    downloadPhaserScript(PHASER_URLS_ONE, '4MRNA-G4-OneBlock.sh', dlPhaserOneBtn)
  );
  dlPhaserTwoBtn?.addEventListener('click', () =>
    downloadPhaserScript(PHASER_URLS_TWO, '4MRNA-G4-TwoBlock.sh', dlPhaserTwoBtn)
  );
});

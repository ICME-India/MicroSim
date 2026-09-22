/**
 * Phase Field 2D/3D Viewer and Animation Controller
 */

const ViewerManager = {
  runId: null,
  metadata: null,
  timesteps: [],
  currentFrameIdx: 0,
  activeField: 'c0',
  activeColormap: 'viridis',
  overlaySlabs: true,
  slabsCount: 4,
  plane3D: 'xy',
  slice3DPos: 0,
  lineCutAxis: 'x',
  lineCutCoord: 64,
  isPlaying: false,
  playTimer: null,
  isLooping: true,
  fps: 10,
  hoverDebounceTimer: null,

  init() {
    this.setupUIListeners();
  },

  setupUIListeners() {
    // Play / Pause
    const btnPlay = document.getElementById('btnPlayPause');
    if (btnPlay) {
      btnPlay.addEventListener('click', () => this.togglePlay());
    }

    // Prev / Next
    const btnPrev = document.getElementById('btnStepBack');
    if (btnPrev) {
      btnPrev.addEventListener('click', () => {
        this.pause();
        this.stepTo(this.currentFrameIdx - 1);
      });
    }
    const btnNext = document.getElementById('btnStepForward');
    if (btnNext) {
      btnNext.addEventListener('click', () => {
        this.pause();
        this.stepTo(this.currentFrameIdx + 1);
      });
    }

    // Loop
    const btnLoop = document.getElementById('btnLoop');
    if (btnLoop) {
      btnLoop.addEventListener('click', () => {
        this.isLooping = !this.isLooping;
        btnLoop.textContent = `🔁 Loop: ${this.isLooping ? 'On' : 'Off'}`;
        btnLoop.classList.toggle('active', this.isLooping);
      });
    }

    // FPS
    const selFps = document.getElementById('playbackFps');
    if (selFps) {
      selFps.addEventListener('change', (e) => {
        this.fps = parseInt(e.target.value) || 10;
        if (this.isPlaying) {
          this.pause();
          this.play();
        }
      });
    }

    // Scrubber
    const scrubber = document.getElementById('timeScrubber');
    if (scrubber) {
      scrubber.addEventListener('input', (e) => {
        this.pause();
        this.stepTo(parseInt(e.target.value));
      });
    }

    // Colormap
    const colormapSel = document.getElementById('colormapSelector');
    if (colormapSel) {
      colormapSel.addEventListener('change', (e) => {
        this.activeColormap = e.target.value;
        this.renderFrame();
      });
    }

    // Slab overlay checkbox & count
    const slabCheck = document.getElementById('overlaySlabsCheck');
    if (slabCheck) {
      slabCheck.addEventListener('change', (e) => {
        this.overlaySlabs = e.target.checked;
        this.renderFrame();
      });
    }
    const slabCount = document.getElementById('overlaySlabsCount');
    if (slabCount) {
      slabCount.addEventListener('change', (e) => {
        this.slabsCount = parseInt(e.target.value) || 2;
        if (this.overlaySlabs) this.renderFrame();
      });
    }

    // 3D Plane & Slice
    const planeSel = document.getElementById('plane3DSelector');
    if (planeSel) {
      planeSel.addEventListener('change', (e) => {
        this.plane3D = e.target.value;
        this.update3DSliceBounds();
        this.renderFrame();
      });
    }
    const sliceSlider = document.getElementById('slice3DSlider');
    if (sliceSlider) {
      sliceSlider.addEventListener('input', (e) => {
        this.slice3DPos = parseInt(e.target.value);
        document.getElementById('slice3DPos').textContent = this.slice3DPos;
        this.renderFrame();
      });
    }

    // Line cut controls
    const btnCutX = document.getElementById('btnCutX');
    const btnCutY = document.getElementById('btnCutY');
    const cutSlider = document.getElementById('lineCutSlider');

    if (btnCutX && btnCutY) {
      btnCutX.addEventListener('click', () => {
        this.lineCutAxis = 'x';
        btnCutX.classList.add('active');
        btnCutY.classList.remove('active');
        document.getElementById('lineCutSliderLabel').textContent = 'Fixed Y:';
        this.updateLineCutBounds();
        this.refreshLineCut();
      });

      btnCutY.addEventListener('click', () => {
        this.lineCutAxis = 'y';
        btnCutY.classList.add('active');
        btnCutX.classList.remove('active');
        document.getElementById('lineCutSliderLabel').textContent = 'Fixed X:';
        this.updateLineCutBounds();
        this.refreshLineCut();
      });
    }

    if (cutSlider) {
      cutSlider.addEventListener('input', (e) => {
        this.lineCutCoord = parseInt(e.target.value);
        document.getElementById('lineCutPos').textContent = this.lineCutCoord;
        document.getElementById('lineCutBadge').textContent = `Cut: ${this.lineCutAxis.toUpperCase()} @ ${this.lineCutCoord}`;
        this.refreshLineCut();
      });
    }

    // Viewer mouse interaction (hover coordinate and click to cut)
    const viewerImg = document.getElementById('viewerImage');
    if (viewerImg) {
      viewerImg.addEventListener('mousemove', (e) => this.handleViewerHover(e));
      viewerImg.addEventListener('mouseleave', () => {
        const hud = document.getElementById('viewerHud');
        if (hud) hud.style.display = 'none';
      });
      viewerImg.addEventListener('click', (e) => this.handleViewerClick(e));
    }

    // Download buttons
    const btnVtk = document.getElementById('btnDownloadVtk');
    if (btnVtk) {
      btnVtk.addEventListener('click', () => {
        if (!this.runId || !this.timesteps[this.currentFrameIdx]) return;
        const step = this.timesteps[this.currentFrameIdx].step;
        window.location.href = `/api/download/vtk/${this.runId}/${step}`;
      });
    }

    const btnPng = document.getElementById('btnDownloadPng');
    if (btnPng) {
      btnPng.addEventListener('click', () => {
        const img = document.getElementById('viewerImage');
        if (img && img.src) {
          const a = document.createElement('a');
          a.href = img.src;
          a.download = `${this.runId}_step_${this.timesteps[this.currentFrameIdx]?.step || 0}.png`;
          document.body.appendChild(a);
          a.click();
          document.body.removeChild(a);
        }
      });
    }

    const btnCsv = document.getElementById('btnDownloadCsv');
    if (btnCsv) {
      btnCsv.addEventListener('click', () => {
        if (!this.runId) return;
        window.location.href = `/api/download/csv/${this.runId}`;
      });
    }
  },

  loadRun(runId, metadata) {
    this.runId = runId;
    this.metadata = metadata;
    this.timesteps = metadata.timesteps || [];
    this.currentFrameIdx = 0;

    // Scrubber bounds
    const scrubber = document.getElementById('timeScrubber');
    if (scrubber) {
      scrubber.min = 0;
      scrubber.max = Math.max(0, this.timesteps.length - 1);
      scrubber.value = 0;
    }

    // Setup component button group
    this.updateComponentButtons(metadata);

    // 3D Controls visibility
    const ctrl3d = document.getElementById('controls3D');
    if (ctrl3d) {
      if (metadata.dim === 3 && metadata.nz > 1) {
        ctrl3d.classList.remove('d-none');
        this.update3DSliceBounds();
      } else {
        ctrl3d.classList.add('d-none');
      }
    }

    // Setup line-cut slider bounds
    this.updateLineCutBounds();

    // Render first frame
    this.stepTo(0);
  },

  updateComponentButtons(metadata) {
    const group = document.getElementById('componentButtonGroup');
    if (!group) return;
    group.innerHTML = '';

    const comps = metadata.components || ['c0', 'c1'];
    comps.forEach((c, idx) => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = `btn btn-ctrl btn-sm ${idx === 0 ? 'active' : ''}`;
      btn.textContent = c;
      btn.addEventListener('click', () => {
        group.querySelectorAll('.btn-ctrl').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        this.activeField = c;
        this.renderFrame();
      });
      group.appendChild(btn);
    });

    // If 3 or more components, offer RGB Composite
    if (comps.length >= 3) {
      const rgbBtn = document.createElement('button');
      rgbBtn.type = 'button';
      rgbBtn.className = 'btn btn-ctrl btn-sm text-info';
      rgbBtn.textContent = 'RGB Composite';
      rgbBtn.title = 'R=c0, G=c1, B=c2 Phase Composite Map';
      rgbBtn.addEventListener('click', () => {
        group.querySelectorAll('.btn-ctrl').forEach(b => b.classList.remove('active'));
        rgbBtn.classList.add('active');
        this.activeField = 'rgb';
        this.renderFrame();
      });
      group.appendChild(rgbBtn);
    }

    // Chemical potentials buttons
    const mus = metadata.chemical_potentials || [];
    mus.forEach(mu => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'btn btn-ctrl btn-sm text-warning';
      btn.textContent = mu;
      btn.addEventListener('click', () => {
        group.querySelectorAll('.btn-ctrl').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        this.activeField = mu;
        this.renderFrame();
      });
      group.appendChild(btn);
    });

    this.activeField = comps[0] || 'c0';
  },

  update3DSliceBounds() {
    if (!this.metadata) return;
    const slider = document.getElementById('slice3DSlider');
    const label = document.getElementById('slice3DLabel');
    if (!slider) return;

    let maxVal = 0;
    if (this.plane3D === 'xy') {
      maxVal = (this.metadata.nz || 1) - 1;
      label.textContent = 'Z Slice:';
    } else if (this.plane3D === 'xz') {
      maxVal = (this.metadata.ny || 1) - 1;
      label.textContent = 'Y Slice:';
    } else if (this.plane3D === 'yz') {
      maxVal = (this.metadata.nx || 1) - 1;
      label.textContent = 'X Slice:';
    }

    slider.min = 0;
    slider.max = Math.max(0, maxVal);
    this.slice3DPos = Math.floor(maxVal / 2);
    slider.value = this.slice3DPos;
    document.getElementById('slice3DPos').textContent = this.slice3DPos;
  },

  updateLineCutBounds() {
    if (!this.metadata) return;
    const slider = document.getElementById('lineCutSlider');
    if (!slider) return;

    if (this.lineCutAxis === 'x') {
      // Cut along X at fixed Y -> max coord is ny-1
      slider.max = (this.metadata.ny || 128) - 1;
    } else {
      // Cut along Y at fixed X -> max coord is nx-1
      slider.max = (this.metadata.nx || 128) - 1;
    }
    this.lineCutCoord = Math.min(this.lineCutCoord, parseInt(slider.max));
    slider.value = this.lineCutCoord;
    document.getElementById('lineCutPos').textContent = this.lineCutCoord;
  },

  stepTo(frameIdx) {
    if (this.timesteps.length === 0) return;
    if (frameIdx < 0) frameIdx = this.isLooping ? this.timesteps.length - 1 : 0;
    if (frameIdx >= this.timesteps.length) frameIdx = this.isLooping ? 0 : this.timesteps.length - 1;

    this.currentFrameIdx = frameIdx;
    const stepInfo = this.timesteps[frameIdx];

    // Update UI counters
    document.getElementById('timeScrubber').value = frameIdx;
    document.getElementById('frameCounter').textContent = `Frame ${frameIdx + 1} / ${this.timesteps.length}`;
    document.getElementById('stepCounter').textContent = `Step: ${stepInfo.step}`;

    this.renderFrame();
    this.refreshLineCut();
  },

  renderFrame() {
    if (!this.runId || this.timesteps.length === 0) return;
    const stepInfo = this.timesteps[this.currentFrameIdx];
    if (!stepInfo) return;

    const img = document.getElementById('viewerImage');
    if (!img) return;

    const slabs = this.overlaySlabs ? this.slabsCount : 0;
    const sliceIdx = (this.metadata?.dim === 3) ? this.slice3DPos : 0;

    const url = `/api/run/${this.runId}/slice_img?step=${stepInfo.step}&field=${this.activeField}&plane=${this.plane3D}&slice_idx=${sliceIdx}&colormap=${this.activeColormap}&overlay_slabs=${slabs}`;

    // Use Fetch to retrieve headers with min/max scale values
    fetch(url)
      .then(res => {
        const vmin = res.headers.get('X-Field-Min') || '0.0000';
        const vmax = res.headers.get('X-Field-Max') || '1.0000';
        const simTime = res.headers.get('X-Time') || '0.000';

        document.getElementById('scaleMin').textContent = `Min: ${vmin}`;
        document.getElementById('scaleMax').textContent = `Max: ${vmax}`;
        document.getElementById('scaleField').textContent = `Field: ${this.activeField}`;
        document.getElementById('timeCounter').textContent = `Time t = ${parseFloat(simTime).toFixed(4)}`;

        return res.blob();
      })
      .then(blob => {
        img.src = URL.createObjectURL(blob);
      })
      .catch(err => console.error('Error rendering slice:', err));
  },

  refreshLineCut() {
    if (!this.runId || this.timesteps.length === 0) return;
    const step = this.timesteps[this.currentFrameIdx]?.step || 0;
    const zSlice = (this.metadata?.dim === 3) ? this.slice3DPos : 0;

    fetch(`/api/run/${this.runId}/linecut?step=${step}&axis=${this.lineCutAxis}&coord=${this.lineCutCoord}&z_slice=${zSlice}`)
      .then(res => res.json())
      .then(data => {
        ChartsManager.updateLineCut(data);
      })
      .catch(err => console.error('Error fetching line cut:', err));
  },

  handleViewerHover(e) {
    const hud = document.getElementById('viewerHud');
    if (!hud || !this.metadata) return;

    const img = e.target;
    const rect = img.getBoundingClientRect();
    const nx = this.metadata.nx || 128;
    const ny = this.metadata.ny || 128;

    const px = Math.floor(((e.clientX - rect.left) / rect.width) * nx);
    const py = Math.floor((1.0 - (e.clientY - rect.top) / rect.height) * ny);

    if (px < 0 || px >= nx || py < 0 || py >= ny) return;

    hud.style.display = 'block';
    document.getElementById('hudCoords').textContent = `Domain Coordinates: (X: ${px}, Y: ${py})`;

    // Debounce hover info fetch
    clearTimeout(this.hoverDebounceTimer);
    this.hoverDebounceTimer = setTimeout(() => {
      const step = this.timesteps[this.currentFrameIdx]?.step || 0;
      fetch(`/api/run/${this.runId}/hover_info?step=${step}&x=${px}&y=${py}`)
        .then(res => res.json())
        .then(data => {
          if (data && data.values) {
            let valStr = Object.entries(data.values)
              .map(([k, v]) => `${k}=${v}`)
              .join(' | ');
            valStr += ` | ∑c=${data.sum_c}`;
            document.getElementById('hudValues').textContent = valStr;
          }
        })
        .catch(() => {});
    }, 50);
  },

  handleViewerClick(e) {
    if (!this.metadata) return;
    const img = e.target;
    const rect = img.getBoundingClientRect();
    const nx = this.metadata.nx || 128;
    const ny = this.metadata.ny || 128;

    const px = Math.floor(((e.clientX - rect.left) / rect.width) * nx);
    const py = Math.floor((1.0 - (e.clientY - rect.top) / rect.height) * ny);

    if (this.lineCutAxis === 'x') {
      this.lineCutCoord = py;
    } else {
      this.lineCutCoord = px;
    }

    const slider = document.getElementById('lineCutSlider');
    if (slider) slider.value = this.lineCutCoord;
    document.getElementById('lineCutPos').textContent = this.lineCutCoord;
    document.getElementById('lineCutBadge').textContent = `Cut: ${this.lineCutAxis.toUpperCase()} @ ${this.lineCutCoord}`;
    this.refreshLineCut();
  },

  play() {
    this.isPlaying = true;
    const btn = document.getElementById('btnPlayPause');
    if (btn) btn.textContent = '⏸ Pause';

    const interval = Math.max(20, Math.floor(1000 / this.fps));
    this.playTimer = setInterval(() => {
      this.stepTo(this.currentFrameIdx + 1);
    }, interval);
  },

  pause() {
    this.isPlaying = false;
    const btn = document.getElementById('btnPlayPause');
    if (btn) btn.textContent = '▶ Play';
    if (this.playTimer) {
      clearInterval(this.playTimer);
      this.playTimer = null;
    }
  },

  togglePlay() {
    if (this.isPlaying) {
      this.pause();
    } else {
      this.play();
    }
  }
};

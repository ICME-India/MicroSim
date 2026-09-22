/**
 * Main Application Controller for Cahn-Hilliard Studio
 */

const App = {
  currentRunId: null,
  pollTimer: null,
  logIndex: 0,
  presets: [],

  init() {
    ChartsManager.initCharts();
    ViewerManager.init();

    this.setupEventListeners();
    this.loadPresets();
    this.loadAvailableRuns();

    // Check status in case simulation was running
    this.pollSimulationStatus();
  },

  setupEventListeners() {
    // Global run selector
    const runSel = document.getElementById('globalRunSelector');
    if (runSel) {
      runSel.addEventListener('change', (e) => {
        this.selectRun(e.target.value);
      });
    }

    // Theme toggle
    const themeBtn = document.getElementById('themeToggleBtn');
    if (themeBtn) {
      themeBtn.addEventListener('click', () => {
        const html = document.documentElement;
        const current = html.getAttribute('data-theme') || 'dark';
        const next = current === 'dark' ? 'light' : 'dark';
        html.setAttribute('data-theme', next);
      });
    }

    // Start simulation button
    const btnStart = document.getElementById('btnStartSim');
    if (btnStart) {
      btnStart.addEventListener('click', () => this.startSimulation());
    }

    // Stop simulation button
    const btnStop = document.getElementById('btnStopSim');
    if (btnStop) {
      btnStop.addEventListener('click', () => this.stopSimulation());
    }

    // Clear logs button
    const btnClear = document.getElementById('btnClearLogs');
    if (btnClear) {
      btnClear.addEventListener('click', () => {
        document.getElementById('consoleLogBox').textContent = '';
      });
    }

    // Tab switch listener to resize charts
    const tabBtns = document.querySelectorAll('#appTabs button');
    tabBtns.forEach(btn => {
      btn.addEventListener('shown.bs.tab', (e) => {
        const targetId = e.target.getAttribute('data-bs-target');
        if (targetId === '#tab-diagnostics') {
          if (this.currentRunId) this.loadDiagnostics(this.currentRunId);
        } else if (targetId === '#tab-visualizer') {
          ViewerManager.renderFrame();
          ViewerManager.refreshLineCut();
        }
      });
    });
  },

  loadAvailableRuns() {
    fetch('/api/runs')
      .then(res => res.json())
      .then(data => {
        const sel = document.getElementById('globalRunSelector');
        if (!sel) return;
        sel.innerHTML = '';

        const runs = data.runs || [];
        if (runs.length === 0) {
          sel.innerHTML = '<option value="">No runs found</option>';
          return;
        }

        runs.forEach(r => {
          const opt = document.createElement('option');
          opt.value = r.id;
          opt.textContent = `${r.title} (${r.vtk_count} frames)`;
          sel.appendChild(opt);
        });

        // Pick preferred initial run
        let initialRun = runs[0].id;
        if (runs.some(r => r.id === 'results_ternary_spinodal')) {
          initialRun = 'results_ternary_spinodal';
        }

        sel.value = initialRun;
        this.selectRun(initialRun);
      })
      .catch(err => console.error('Error loading runs:', err));
  },

  selectRun(runId) {
    if (!runId) return;
    this.currentRunId = runId;

    // Fetch metadata
    fetch(`/api/run/${runId}/metadata`)
      .then(res => res.json())
      .then(meta => {
        // Update metadata card
        document.getElementById('metaGrid').textContent = `${meta.nx} × ${meta.ny}` + (meta.dim === 3 ? ` × ${meta.nz}` : '');
        document.getElementById('metaComponents').textContent = (meta.components || []).join(', ');
        document.getElementById('metaModel').textContent = meta.config?.free_energy_type || 'polynomial_multiwell';
        document.getElementById('metaIntegrator').textContent = meta.config?.integrator || 'rk2';
        document.getElementById('metaDt').textContent = meta.config?.dt || '0.01';
        document.getElementById('metaFrames').textContent = meta.total_frames || '0';

        // Update Viewer
        ViewerManager.loadRun(runId, meta);

        // Update Diagnostics
        this.loadDiagnostics(runId);
      })
      .catch(err => console.error('Error fetching run metadata:', err));
  },

  loadDiagnostics(runId) {
    fetch(`/api/run/${runId}/diagnostics`)
      .then(res => {
        if (!res.ok) throw new Error('No diagnostics available');
        return res.json();
      })
      .then(diag => {
        if (!diag || diag.error) return;

        // KPI cards
        document.getElementById('diagInitialF').textContent = diag.initial_free_energy.toFixed(3);
        document.getElementById('diagFinalF').textContent = `Final: ${diag.final_free_energy.toFixed(3)}`;
        document.getElementById('diagPctDrop').textContent = `${Math.abs(diag.pct_free_energy_decrease).toFixed(2)}%`;
        document.getElementById('diagMaxDrift').textContent = diag.max_mass_drift.toExponential(3);
        document.getElementById('diagUnityDev').textContent = diag.max_unity_deviation.toExponential(3);

        const checkElem = document.getElementById('diagDissipCheck');
        if (diag.is_monotone_dissipative) {
          checkElem.textContent = 'dF/dt ≤ 0 Verified (Monotonic)';
          checkElem.className = 'kpi-sub text-success';
        } else {
          checkElem.textContent = 'dF/dt check warning';
          checkElem.className = 'kpi-sub text-warning';
        }

        ChartsManager.updateDiagnostics(diag);
      })
      .catch(err => {
        console.log('Diagnostics not loaded:', err);
      });
  },

  loadPresets() {
    fetch('/api/presets')
      .then(res => res.json())
      .then(data => {
        this.presets = data.presets || [];
        const container = document.getElementById('presetCardsContainer');
        if (!container) return;
        container.innerHTML = '';

        this.presets.forEach((p, idx) => {
          const col = document.createElement('div');
          col.className = 'col-md-4 col-lg-2';
          col.innerHTML = `
            <div class="preset-card ${idx === 1 ? 'active' : ''}" data-file="${p.file}">
              <div class="fw-bold" style="font-size: 0.88rem; color: #38bdf8;">${p.title}</div>
              <div class="text-secondary" style="font-size: 0.75rem; margin-top: 4px;">${p.description}</div>
            </div>
          `;
          col.querySelector('.preset-card').addEventListener('click', (e) => {
            document.querySelectorAll('.preset-card').forEach(c => c.classList.remove('active'));
            col.querySelector('.preset-card').classList.add('active');
            this.populateConfigForm(p.config, p.file);
          });
          container.appendChild(col);
        });

        // Initialize with ternary spinodal preset if available
        const defaultPreset = this.presets.find(p => p.file === 'ternary_spinodal.json') || this.presets[0];
        if (defaultPreset) {
          this.populateConfigForm(defaultPreset.config, defaultPreset.file);
        }
      });
  },

  populateConfigForm(cfg, filename) {
    if (!cfg) return;

    document.getElementById('cfgDim').value = cfg.dim || 2;
    document.getElementById('cfgNx').value = cfg.nx || 128;
    document.getElementById('cfgNy').value = cfg.ny || 128;
    document.getElementById('cfgNz').value = cfg.nz || 1;
    document.getElementById('cfgDx').value = cfg.dx || 1.0;
    document.getElementById('cfgDy').value = cfg.dy || 1.0;
    document.getElementById('cfgBcX').value = cfg.bc_x || 'periodic';
    document.getElementById('cfgBcY').value = cfg.bc_y || 'periodic';

    document.getElementById('cfgFreeEnergy').value = cfg.free_energy_type || 'polynomial_multiwell';
    document.getElementById('cfgNumComp').value = cfg.num_components || 3;
    document.getElementById('cfgMobilityType').value = cfg.mobility_type || 'constant';
    let mobVal = cfg.mobility_val;
    if (cfg.mobility && Array.isArray(cfg.mobility) && cfg.mobility.length > 0) {
      mobVal = Array.isArray(cfg.mobility[0]) ? cfg.mobility[0][0] : cfg.mobility[0];
    } else if (cfg.mobility_matrix && Array.isArray(cfg.mobility_matrix) && cfg.mobility_matrix.length > 0) {
      mobVal = Array.isArray(cfg.mobility_matrix[0]) ? cfg.mobility_matrix[0][0] : cfg.mobility_matrix[0];
    }
    document.getElementById('cfgMobilityVal').value = mobVal !== undefined ? mobVal : 1.0;

    document.getElementById('cfgIntegrator').value = cfg.integrator || 'rk2';
    document.getElementById('cfgDt').value = cfg.dt || 0.01;
    document.getElementById('cfgTotalSteps').value = cfg.total_steps || 1500;
    document.getElementById('cfgOutputInterval').value = cfg.output_interval || 150;
    document.getElementById('cfgDiagInterval').value = cfg.diag_interval || 50;

    document.getElementById('cfgInitialCondition').value = cfg.initial_condition || 'random';
    if (cfg.c_mean) {
      document.getElementById('cfgCMean').value = cfg.c_mean.join(', ');
    }
    document.getElementById('cfgNoiseAmp').value = cfg.noise_amp !== undefined ? cfg.noise_amp : 0.05;
    document.getElementById('cfgSeed').value = cfg.seed || 42;

    const baseName = (filename || 'custom_run').replace('.json', '');
    document.getElementById('cfgOutputDir').value = `results_${baseName}`;
  },

  buildConfigFromForm() {
    const dim = parseInt(document.getElementById('cfgDim').value);
    const numComp = parseInt(document.getElementById('cfgNumComp').value);
    const feType = document.getElementById('cfgFreeEnergy').value;
    const cMeanStr = document.getElementById('cfgCMean').value;
    const cMean = cMeanStr.split(',').map(s => parseFloat(s.trim())).filter(n => !isNaN(n));

    const config = {
      dim: dim,
      nx: parseInt(document.getElementById('cfgNx').value),
      ny: parseInt(document.getElementById('cfgNy').value),
      nz: parseInt(document.getElementById('cfgNz').value),
      dx: parseFloat(document.getElementById('cfgDx').value),
      dy: parseFloat(document.getElementById('cfgDy').value),
      bc_x: document.getElementById('cfgBcX').value,
      bc_y: document.getElementById('cfgBcY').value,
      num_components: numComp,
      free_energy_type: feType,
      mobility_type: document.getElementById('cfgMobilityType').value,
      mobility_val: parseFloat(document.getElementById('cfgMobilityVal').value),
      integrator: document.getElementById('cfgIntegrator').value,
      dt: parseFloat(document.getElementById('cfgDt').value),
      total_steps: parseInt(document.getElementById('cfgTotalSteps').value),
      output_interval: parseInt(document.getElementById('cfgOutputInterval').value),
      diag_interval: parseInt(document.getElementById('cfgDiagInterval').value),
      initial_condition: document.getElementById('cfgInitialCondition').value,
      c_mean: cMean.length === numComp ? cMean : Array(numComp).fill(1.0 / numComp),
      noise_amp: parseFloat(document.getElementById('cfgNoiseAmp').value),
      seed: parseInt(document.getElementById('cfgSeed').value),
      output_dir: document.getElementById('cfgOutputDir').value || 'results_custom'
    };

    if (numComp > 2 && config.mobility_type === 'constant') {
      const m0 = config.mobility_val;
      config.mobility = Array.from({length: numComp}, (_, i) =>
        Array.from({length: numComp}, (_, j) => (i === j ? m0 : 0.0))
      );
    }

    if (config.initial_condition === 'droplet') {
      config.droplet_comp = 0;
      config.droplet_radius = 28.0;
      config.droplet_diffuse_width = 2.5;
    }

    return config;
  },

  startSimulation() {
    const config = this.buildConfigFromForm();
    const backend = document.getElementById('cfgBackend').value;
    const device = document.getElementById('cfgDevice') ? document.getElementById('cfgDevice').value : 'cpu';
    const mpiRanks = parseInt(document.getElementById('cfgMpiRanks').value) || 2;

    fetch('/api/simulation/start', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        config: config,
        backend: backend,
        device: device,
        mpi_ranks: mpiRanks
      })
    })
    .then(res => res.json())
    .then(data => {
      if (data.success) {
        this.logIndex = 0;
        const devLabel = device.toUpperCase() === 'GPU' ? 'GPU (CUDA)' : 'CPU';
        document.getElementById('consoleLogBox').textContent = `[STUDIO] Starting ${backend.toUpperCase()} simulation on ${devLabel} with ${mpiRanks} ranks...\n`;
        this.pollSimulationStatus();
      } else {
        alert('Failed to start simulation: ' + (data.error || 'Unknown error'));
      }
    })
    .catch(err => alert('Network error: ' + err));
  },

  stopSimulation() {
    fetch('/api/simulation/stop', { method: 'POST' })
      .then(res => res.json())
      .then(data => {
        console.log('Stop response:', data);
      });
  },

  pollSimulationStatus() {
    clearInterval(this.pollTimer);

    const poll = () => {
      fetch('/api/simulation/status')
        .then(res => res.json())
        .then(status => {
          this.updateStatusUI(status);

          // Fetch new logs
          fetch(`/api/simulation/logs?since=${this.logIndex}`)
            .then(r => r.json())
            .then(logData => {
              if (logData.lines && logData.lines.length > 0) {
                const box = document.getElementById('consoleLogBox');
                box.textContent += logData.lines.join('\n') + '\n';
                box.scrollTop = box.scrollHeight;
                this.logIndex = logData.next_index;
              }
            });

          if (status.status === 'running') {
            // Keep polling
          } else {
            // Finished, stopped, or idle
            clearInterval(this.pollTimer);
            if (status.status === 'completed') {
              this.loadAvailableRuns();
              // Switch to newly generated run
              if (status.run_id) {
                setTimeout(() => {
                  const sel = document.getElementById('globalRunSelector');
                  if (sel) {
                    sel.value = status.run_id;
                    this.selectRun(status.run_id);
                  }
                }, 800);
              }
            }
          }
        })
        .catch(() => clearInterval(this.pollTimer));
    };

    poll();
    this.pollTimer = setInterval(poll, 600);
  },

  updateStatusUI(status) {
    const pill = document.getElementById('statusPill');
    const pillText = document.getElementById('statusPillText');
    const badge = document.getElementById('simStateBadge');
    const btnStart = document.getElementById('btnStartSim');
    const btnStop = document.getElementById('btnStopSim');
    const progressBar = document.getElementById('simProgressBar');

    pill.className = `status-pill status-${status.status}`;
    pillText.textContent = status.status.toUpperCase();
    badge.className = `badge bg-${status.status === 'running' ? 'success' : (status.status === 'completed' ? 'primary' : 'secondary')}`;
    badge.textContent = status.status.toUpperCase();

    if (status.status === 'running') {
      btnStart.disabled = true;
      btnStop.disabled = false;
    } else {
      btnStart.disabled = false;
      btnStop.disabled = true;
    }

    progressBar.style.width = `${status.progress_percent}%`;
    progressBar.textContent = `${status.progress_percent}%`;

    document.getElementById('liveStep').textContent = `${status.current_step} / ${status.total_steps}`;
    document.getElementById('liveTime').textContent = status.current_time ? status.current_time.toFixed(4) : '0.0000';
    document.getElementById('liveEnergy').textContent = status.current_energy ? status.current_energy.toExponential(3) : '--';
    document.getElementById('liveMlups').textContent = status.current_mlups ? `${status.current_mlups.toFixed(2)} MLUPS` : '-- MLUPS';
  }
};

// Start application on DOM loaded
document.addEventListener('DOMContentLoaded', () => {
  App.init();
});

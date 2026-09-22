/**
 * Interactive Diagnostic and Line-Cut Charts using Chart.js
 */

const ChartsManager = {
  freeEnergyChart: null,
  dfDtChart: null,
  massDriftChart: null,
  unityDevChart: null,
  lineCutChart: null,

  chartDefaults() {
    return {
      responsive: true,
      maintainAspectRatio: false,
      animation: { duration: 250 },
      plugins: {
        legend: {
          labels: { color: '#9ca3af', font: { size: 11 } }
        },
        tooltip: {
          mode: 'index',
          intersect: false
        }
      },
      scales: {
        x: {
          grid: { color: 'rgba(255, 255, 255, 0.05)' },
          ticks: { color: '#9ca3af', maxTicksLimit: 8 }
        },
        y: {
          grid: { color: 'rgba(255, 255, 255, 0.05)' },
          ticks: { color: '#9ca3af', maxTicksLimit: 6 }
        }
      }
    };
  },

  initCharts() {
    const isChartAvailable = typeof Chart !== 'undefined';
    if (!isChartAvailable) {
      console.warn('Chart.js not loaded yet');
      return;
    }

    // Line Cut Chart
    const ctxLine = document.getElementById('lineCutChart');
    if (ctxLine && !this.lineCutChart) {
      this.lineCutChart = new Chart(ctxLine, {
        type: 'line',
        data: { labels: [], datasets: [] },
        options: {
          ...this.chartDefaults(),
          plugins: {
            legend: { display: true, labels: { color: '#cbd5e1' } }
          }
        }
      });
    }

    // Diagnostics Charts
    const ctxF = document.getElementById('chartFreeEnergy');
    if (ctxF && !this.freeEnergyChart) {
      this.freeEnergyChart = new Chart(ctxF, {
        type: 'line',
        data: { labels: [], datasets: [] },
        options: this.chartDefaults()
      });
    }

    const ctxDf = document.getElementById('chartDfDt');
    if (ctxDf && !this.dfDtChart) {
      this.dfDtChart = new Chart(ctxDf, {
        type: 'line',
        data: { labels: [], datasets: [] },
        options: this.chartDefaults()
      });
    }

    const ctxMass = document.getElementById('chartMassDrift');
    if (ctxMass && !this.massDriftChart) {
      this.massDriftChart = new Chart(ctxMass, {
        type: 'line',
        data: { labels: [], datasets: [] },
        options: this.chartDefaults()
      });
    }

    const ctxDev = document.getElementById('chartUnityDev');
    if (ctxDev && !this.unityDevChart) {
      this.unityDevChart = new Chart(ctxDev, {
        type: 'line',
        data: { labels: [], datasets: [] },
        options: {
          ...this.chartDefaults(),
          scales: {
            ...this.chartDefaults().scales,
            y: {
              ...this.chartDefaults().scales.y,
              type: 'logarithmic'
            }
          }
        }
      });
    }
  },

  updateLineCut(data) {
    if (!this.lineCutChart || !data || !data.positions) return;

    const colors = ['#38bdf8', '#34d399', '#f43f5e', '#fbbf24', '#a78bfa', '#f97316'];
    const datasets = [];

    let cIdx = 0;
    for (const [key, profile] of Object.entries(data.profiles || {})) {
      datasets.push({
        label: key,
        data: profile,
        borderColor: colors[cIdx % colors.length],
        backgroundColor: colors[cIdx % colors.length] + '22',
        borderWidth: 2,
        pointRadius: 0,
        fill: false,
        tension: 0.1
      });
      cIdx++;
    }

    this.lineCutChart.data.labels = data.positions;
    this.lineCutChart.data.datasets = datasets;
    this.lineCutChart.options.scales.x.title = {
      display: true,
      text: data.axis === 'x' ? 'X Position' : 'Y Position',
      color: '#9ca3af'
    };
    this.lineCutChart.update();
  },

  updateDiagnostics(diag) {
    if (!diag || !diag.series) return;

    const series = diag.series;
    const times = series.times;

    // 1. Free energy chart
    if (this.freeEnergyChart) {
      this.freeEnergyChart.data.labels = times;
      this.freeEnergyChart.data.datasets = [{
        label: 'Total Free Energy F(t)',
        data: series.total_free_energy,
        borderColor: '#38bdf8',
        backgroundColor: 'rgba(56, 189, 248, 0.1)',
        borderWidth: 2,
        pointRadius: 0,
        fill: true,
        tension: 0.1
      }];
      this.freeEnergyChart.update();
    }

    // 2. dF/dt chart
    if (this.dfDtChart) {
      this.dfDtChart.data.labels = times;
      this.dfDtChart.data.datasets = [{
        label: 'Rate dF/dt',
        data: series.dF_dt,
        borderColor: '#f43f5e',
        backgroundColor: 'rgba(244, 63, 94, 0.1)',
        borderWidth: 1.8,
        pointRadius: 0,
        tension: 0.1
      }];
      this.dfDtChart.update();
    }

    // 3. Mass drift chart
    if (this.massDriftChart) {
      const colors = ['#34d399', '#fbbf24', '#a78bfa', '#f97316', '#38bdf8'];
      const datasets = [];
      let cIdx = 0;
      for (const [comp, drifts] of Object.entries(series.mass_drifts || {})) {
        datasets.push({
          label: `Δc̄_${comp.replace('avg_c', '')}`,
          data: drifts,
          borderColor: colors[cIdx % colors.length],
          borderWidth: 1.5,
          pointRadius: 0,
          tension: 0.1
        });
        cIdx++;
      }
      this.massDriftChart.data.labels = times;
      this.massDriftChart.data.datasets = datasets;
      this.massDriftChart.update();
    }

    // 4. Partition of unity chart
    if (this.unityDevChart) {
      this.unityDevChart.data.labels = times;
      this.unityDevChart.data.datasets = [{
        label: 'Max |∑ cᵢ - 1|',
        data: series.unity_deviation.map(v => Math.max(1e-18, v)),
        borderColor: '#fbbf24',
        borderWidth: 1.8,
        pointRadius: 0,
        tension: 0.1
      }];
      this.unityDevChart.update();
    }
  }
};

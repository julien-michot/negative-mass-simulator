const canvas = document.getElementById('gpuCanvas');
const statusEl = document.getElementById('status');

let PARTICLE_COUNT = 5000;
const GRID_LEVELS = [32, 16, 8, 4, 2, 1];
const LEVEL_COUNT = GRID_LEVELS.length;
const DOMAIN_SIZE = 1.0;
const DT = 0.00002;
const THETA = 0.5;
const SOFTENING = 0.001;
const DAMPING = 1.0;
const GRAVITY = 10.0;
const MASS_RANGE = [0.5, 1.0];
const MASS_SCALE = 1000; // fixed-point scale for atomic mass
const POS_SCALE = 10000; // fixed-point scale for center-of-mass accumulation (increase to reduce quantization)
const USE_MERGE_GRID = false;
const MERGE_GRID = USE_MERGE_GRID ? 128 : 1;
const MERGE_REACH = 0.012;
const MERGE_POS_SCALE = 200; // lower scale for merge accumulation to avoid u32 overflow
const MERGE_VEL_SCALE = 1000; // lower scale for merge velocity accumulation to avoid i32 overflow
const VEL_SCALE = 10000; // increase velocity scale to reduce quantization in merge accumulators
const ENABLE_MERGING = false; // overall merge toggle; set false to avoid aggressive merging
const USE_DIRECT_SUM = true; // exact O(N^2) forces for smooth, noise-free motion

const uiState = {
  dt: DT,
  gravity: GRAVITY,
  theta: THETA,
  useDirect: USE_DIRECT_SUM
};
// visibility flags for positive / negative particles
uiState.showPositive = true;
uiState.showNegative = true;
// runtime merge toggle (UI-controlled)
uiState.enableMerging = ENABLE_MERGING;
// runtime merge reach (UI-controlled)
uiState.mergeReach = MERGE_REACH;

// runtime pause flag (toggled by UI)
let isPaused = false;
// request a single compute step while paused
let stepRequested = false;
// if user changes solver before regenerateParticles is ready, defer reset
let resetOnMethodChange = false;

const WORKGROUP_SIZE = 8;
const WORKGROUP_SIZE_1D = 256;

const viewState = {
  center: [0.5, 0.5],
  zoom: 0.5
};

const pointerState = {
  active: new Map(),
  lastPan: null,
  lastPinchDist: null
};

function setStatus(message) {
  statusEl.textContent = message;
}

function setupControls() {
  const solverMode = document.getElementById('solverMode');
  const dtSlider = document.getElementById('dtSlider');
  const gravitySlider = document.getElementById('gravitySlider');
  const thetaSlider = document.getElementById('thetaSlider');
  const particleCountSlider = document.getElementById('particleCountSlider');
  const particleCountValue = document.getElementById('particleCountValue');
  const regenerateButton = document.getElementById('regenerateButton');
  const dtValue = document.getElementById('dtValue');
  const gravityValue = document.getElementById('gravityValue');
  const thetaValue = document.getElementById('thetaValue');
  const showPositiveCheckbox = document.getElementById('showPositive');
  const showNegativeCheckbox = document.getElementById('showNegative');
  const enableMergingCheckbox = document.getElementById('enableMerging');

  if (!solverMode || !dtSlider || !gravitySlider || !thetaSlider) {
    return;
  }

  solverMode.value = uiState.useDirect ? 'direct' : 'barnes';
  // dtSlider is log-scaled: slider value represents exponent (base 10)
  // default to log10(uiState.dt) when possible, otherwise fall back to -3
  dtSlider.value = (uiState.dt > 0 ? Math.log10(uiState.dt) : -3).toString();
  gravitySlider.value = uiState.gravity.toString();
  thetaSlider.value = uiState.theta.toString();
  if (particleCountSlider) particleCountSlider.value = String(PARTICLE_COUNT);
  if (particleCountValue) particleCountValue.textContent = String(PARTICLE_COUNT);

  if (showPositiveCheckbox) showPositiveCheckbox.checked = uiState.showPositive;
  if (showNegativeCheckbox) showNegativeCheckbox.checked = uiState.showNegative;
  if (enableMergingCheckbox) enableMergingCheckbox.checked = uiState.enableMerging;

  const updateLabels = () => {
    if (dtValue) dtValue.textContent = Number(uiState.dt).toExponential(3);
    if (gravityValue) gravityValue.textContent = Number(uiState.gravity).toFixed(2);
    if (thetaValue) thetaValue.textContent = Number(uiState.theta).toFixed(2);
  };

  solverMode.addEventListener('change', () => {
    uiState.useDirect = solverMode.value === 'direct';
    // reset the simulation when switching solver to avoid stale buffer state
    if (window.regenerateParticles) {
      window.regenerateParticles(PARTICLE_COUNT);
    } else {
      resetOnMethodChange = true;
    }
  });
  dtSlider.addEventListener('input', () => {
    // slider gives exponent; compute dt = 10^exponent
    const exp = Number(dtSlider.value);
    uiState.dt = Math.pow(10, exp);
    updateLabels();
  });

  const pauseButton = document.getElementById('pauseButton');
  if (pauseButton) {
    pauseButton.textContent = isPaused ? 'Resume' : 'Pause';
    pauseButton.addEventListener('click', () => {
      isPaused = !isPaused;
      pauseButton.textContent = isPaused ? 'Resume' : 'Pause';
      setStatus(isPaused ? 'Paused' : 'Running Barnes-Hut solver on GPU.');
    });
  }
  const stepButton = document.getElementById('stepButton');
  if (stepButton) {
    stepButton.addEventListener('click', () => {
      // request a single compute step; leave isPaused state unchanged
      stepRequested = true;
      setStatus('Single step requested');
    });
  }
  gravitySlider.addEventListener('input', () => {
    uiState.gravity = Number(gravitySlider.value);
    updateLabels();
  });
  thetaSlider.addEventListener('input', () => {
    uiState.theta = Number(thetaSlider.value);
    updateLabels();
  });

  if (showPositiveCheckbox) {
    showPositiveCheckbox.addEventListener('change', () => {
      uiState.showPositive = !!showPositiveCheckbox.checked;
    });
  }
  if (showNegativeCheckbox) {
    showNegativeCheckbox.addEventListener('change', () => {
      uiState.showNegative = !!showNegativeCheckbox.checked;
    });
  }
  if (enableMergingCheckbox) {
    enableMergingCheckbox.addEventListener('change', () => {
      uiState.enableMerging = !!enableMergingCheckbox.checked;
    });
  }

  if (particleCountSlider) {
    particleCountSlider.addEventListener('input', () => {
      if (particleCountValue) particleCountValue.textContent = String(particleCountSlider.value);
    });
  }

  if (regenerateButton && particleCountSlider) {
    regenerateButton.addEventListener('click', () => {
      const count = Number(particleCountSlider.value) || PARTICLE_COUNT;
      if (window.regenerateParticles) {
        window.regenerateParticles(count);
      } else {
        setStatus('Regeneration not ready');
      }
      if (particleCountValue) particleCountValue.textContent = String(count);
    });
  }

  // merge reach slider wiring
  const mergeReachSlider = document.getElementById('mergeReachSlider');
  const mergeReachValue = document.getElementById('mergeReachValue');
  if (mergeReachSlider) {
    mergeReachSlider.value = String(uiState.mergeReach);
    if (mergeReachValue) mergeReachValue.textContent = Number(uiState.mergeReach).toFixed(2);
    mergeReachSlider.addEventListener('input', () => {
      uiState.mergeReach = Number(mergeReachSlider.value);
      if (mergeReachValue) mergeReachValue.textContent = Number(uiState.mergeReach).toFixed(2);
    });
  }

  updateLabels();
}

function createBuffer(device, data, usage) {
  const buffer = device.createBuffer({
    size: data.byteLength,
    usage,
    mappedAtCreation: true
  });
  const constructor = data.constructor;
  new constructor(buffer.getMappedRange()).set(data);
  buffer.unmap();
  return buffer;
}

function createEmptyBuffer(device, size, usage) {
  return device.createBuffer({ size, usage });
}

async function loadTexture(device, url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Failed to load texture: ${response.status}`);
  }
  const blob = await response.blob();
  const bitmap = await createImageBitmap(blob);
  const texture = device.createTexture({
    size: [bitmap.width, bitmap.height, 1],
    format: 'rgba8unorm',
    usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT
  });
  device.queue.copyExternalImageToTexture(
    { source: bitmap },
    { texture },
    [bitmap.width, bitmap.height]
  );
  return texture;
}

function buildSimParams() {
  const params = new Float32Array(16);
  params[0] = LEVEL_COUNT;
  params[1] = PARTICLE_COUNT;
  params[2] = DOMAIN_SIZE;
  params[3] = uiState.dt;
  params[4] = uiState.theta;
  params[5] = SOFTENING;
  params[6] = DAMPING;
  params[7] = uiState.gravity;
  params[8] = MASS_SCALE;
  params[9] = 1 / MASS_SCALE;
  params[10] = POS_SCALE;
  params[11] = 1 / POS_SCALE;
  // Disable merge skipping when merging is off; otherwise use UI-controlled mergeReach
  params[12] = uiState.enableMerging ? uiState.mergeReach : 0.0;
  params[13] = VEL_SCALE;
  params[14] = 1 / VEL_SCALE;
  return params;
}

function buildViewParams(aspect) {
  const params = new Float32Array(8);
  params[0] = viewState.center[0];
  params[1] = viewState.center[1];
  params[2] = viewState.zoom;
  params[3] = aspect;
  // showPositive / showNegative flags (1.0 = show, 0.0 = hide)
  params[4] = uiState.showPositive ? 1.0 : 0.0;
  params[5] = uiState.showNegative ? 1.0 : 0.0;
  return params;
}

function buildLevelInfo() {
  const info = new Uint32Array(LEVEL_COUNT * 2);
  let offset = 0;
  for (let i = 0; i < LEVEL_COUNT; i += 1) {
    info[i * 2] = GRID_LEVELS[i];
    info[i * 2 + 1] = offset;
    offset += GRID_LEVELS[i] * GRID_LEVELS[i];
  }
  return { info, totalCells: offset };
}

function createParticles(count = PARTICLE_COUNT) {
  const data = new Float32Array(count * 6);
  for (let i = 0; i < count; i += 1) {
    const base = i * 6;
  const angle = Math.random() * Math.PI * 2;
  const radius = Math.sqrt(Math.random()) * 0.4;
    const cx = 0.5 + Math.cos(angle) * radius;
    const cy = 0.5 + Math.sin(angle) * radius;
    const massMag = MASS_RANGE[0] + Math.random() * (MASS_RANGE[1] - MASS_RANGE[0]);
    const isNegative = i < Math.floor(count / 2);
    data[base] = cx;
    data[base + 1] = cy;
    data[base + 2] = 0;
    data[base + 3] = 0;
    data[base + 4] = isNegative ? -massMag : massMag;
    data[base + 5] = 0;
  }
  return data;
}

function setupInteractions() {
  canvas.addEventListener('wheel', event => {
    event.preventDefault();
    const zoomFactor = Math.exp(-event.deltaY * 0.0015);
    viewState.zoom = Math.min(Math.max(viewState.zoom * zoomFactor, 0.05), 4.0);
  }, { passive: false });

  canvas.addEventListener('pointerdown', event => {
    canvas.setPointerCapture(event.pointerId);
    pointerState.active.set(event.pointerId, { x: event.clientX, y: event.clientY });
    if (pointerState.active.size === 1) {
      pointerState.lastPan = { x: event.clientX, y: event.clientY };
    } else if (pointerState.active.size === 2) {
      const points = [...pointerState.active.values()];
      pointerState.lastPinchDist = Math.hypot(points[0].x - points[1].x, points[0].y - points[1].y);
    }
  });

  canvas.addEventListener('pointermove', event => {
    if (!pointerState.active.has(event.pointerId)) return;
    pointerState.active.set(event.pointerId, { x: event.clientX, y: event.clientY });

    if (pointerState.active.size === 1 && pointerState.lastPan) {
      const dx = event.clientX - pointerState.lastPan.x;
      const dy = event.clientY - pointerState.lastPan.y;
      pointerState.lastPan = { x: event.clientX, y: event.clientY };
      viewState.center[0] -= dx / canvas.clientWidth / viewState.zoom;
      viewState.center[1] += dy / canvas.clientHeight / viewState.zoom;
    } else if (pointerState.active.size === 2) {
      const points = [...pointerState.active.values()];
      const dist = Math.hypot(points[0].x - points[1].x, points[0].y - points[1].y);
      if (pointerState.lastPinchDist) {
        const zoomFactor = dist / pointerState.lastPinchDist;
        viewState.zoom = Math.min(Math.max(viewState.zoom * zoomFactor, 0.005), 4.0);
      }
      pointerState.lastPinchDist = dist;
    }
  });

  const clearPointer = event => {
    pointerState.active.delete(event.pointerId);
    if (pointerState.active.size === 0) {
      pointerState.lastPan = null;
      pointerState.lastPinchDist = null;
    }
  };

  canvas.addEventListener('pointerup', clearPointer);
  canvas.addEventListener('pointercancel', clearPointer);
}

const shaderSources = {
  clearLevels: `
struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read_write> levelMass: array<atomic<i32>>;
@group(0) @binding(1) var<storage, read_write> levelSumX: array<atomic<i32>>;
@group(0) @binding(2) var<storage, read_write> levelSumY: array<atomic<i32>>;
@group(0) @binding(3) var<uniform> params: SimParams;
@group(0) @binding(4) var<storage, read> levelInfo: array<vec2<u32>>;

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let idx = gid.x;
  let levelCount = u32(params.levelCount);
  let last = levelInfo[(levelCount - 1u)].y + levelInfo[(levelCount - 1u)].x * levelInfo[(levelCount - 1u)].x;
  if (idx >= last) {
    return;
  }
  atomicStore(&levelMass[idx], 0);
  atomicStore(&levelSumX[idx], 0);
  atomicStore(&levelSumY[idx], 0);
}
`,
  deposit: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read_write> levelMass: array<atomic<i32>>;
@group(0) @binding(1) var<storage, read_write> levelSumX: array<atomic<i32>>;
@group(0) @binding(2) var<storage, read_write> levelSumY: array<atomic<i32>>;
@group(0) @binding(3) var<storage, read> particles: array<Particle>;
@group(0) @binding(4) var<uniform> params: SimParams;
@group(0) @binding(5) var<storage, read> levelInfo: array<vec2<u32>>;

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let count = u32(params.particleCount);
  if (gid.x >= count) {
    return;
  }
  let p = particles[gid.x];
  if (p.mass == 0.0) {
    return;
  }
  if (p.pos.x < 0.0 || p.pos.x > 1.0 || p.pos.y < 0.0 || p.pos.y > 1.0) {
    return;
  }
  let size = i32(levelInfo[0].x);
  let fx = p.pos.x * f32(size);
  let fy = p.pos.y * f32(size);
  let x = i32(clamp(floor(fx), 0.0, f32(size - 1)));
  let y = i32(clamp(floor(fy), 0.0, f32(size - 1)));
  let idx = u32(y * size + x);

  let massScaled = i32(round(p.mass * params.massScale));
  let sumX = i32(round(p.pos.x * params.posScale * f32(massScaled)));
  let sumY = i32(round(p.pos.y * params.posScale * f32(massScaled)));

  atomicAdd(&levelMass[idx], massScaled);
  atomicAdd(&levelSumX[idx], sumX);
  atomicAdd(&levelSumY[idx], sumY);
}
`,
  reduceLevels: `
struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read_write> levelMass: array<atomic<i32>>;
@group(0) @binding(1) var<storage, read_write> levelSumX: array<atomic<i32>>;
@group(0) @binding(2) var<storage, read_write> levelSumY: array<atomic<i32>>;
@group(0) @binding(3) var<uniform> params: SimParams;
@group(0) @binding(4) var<storage, read> levelInfo: array<vec2<u32>>;

@compute @workgroup_size(${WORKGROUP_SIZE}, ${WORKGROUP_SIZE})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let level = u32(params.dt);
  if (level == 0u) {
    return;
  }
  let parentSize = levelInfo[level].x;
  if (gid.x >= parentSize || gid.y >= parentSize) {
    return;
  }
  let childSize = levelInfo[level - 1u].x;
  let parentOffset = levelInfo[level].y;
  let childOffset = levelInfo[level - 1u].y;

  let baseX = gid.x * 2u;
  let baseY = gid.y * 2u;

  let idx00 = childOffset + baseY * childSize + baseX;
  let idx10 = childOffset + baseY * childSize + baseX + 1u;
  let idx01 = childOffset + (baseY + 1u) * childSize + baseX;
  let idx11 = childOffset + (baseY + 1u) * childSize + baseX + 1u;

  let mass = atomicLoad(&levelMass[idx00]) + atomicLoad(&levelMass[idx10]) + atomicLoad(&levelMass[idx01]) + atomicLoad(&levelMass[idx11]);
  let sumX = atomicLoad(&levelSumX[idx00]) + atomicLoad(&levelSumX[idx10]) + atomicLoad(&levelSumX[idx01]) + atomicLoad(&levelSumX[idx11]);
  let sumY = atomicLoad(&levelSumY[idx00]) + atomicLoad(&levelSumY[idx10]) + atomicLoad(&levelSumY[idx01]) + atomicLoad(&levelSumY[idx11]);

  let parentIdx = parentOffset + gid.y * parentSize + gid.x;
  atomicStore(&levelMass[parentIdx], mass);
  atomicStore(&levelSumX[parentIdx], sumX);
  atomicStore(&levelSumY[parentIdx], sumY);
}
`,
  updateParticles: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read_write> particles: array<Particle>;
@group(0) @binding(1) var<storage, read> levelInfo: array<vec2<u32>>;
@group(0) @binding(2) var<storage, read_write> levelMass: array<atomic<i32>>;
@group(0) @binding(3) var<storage, read_write> levelSumX: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> levelSumY: array<atomic<i32>>;
@group(0) @binding(5) var<uniform> params: SimParams;

fn wrapDelta(delta: f32) -> f32 {
  // No periodic wrapping: use straight delta for open/infinite space
  return delta;
}

fn cellCenter(mass: i32, sumX: i32, sumY: i32, invPosScale: f32) -> vec2<f32> {
  let massF = f32(mass);
  let center = vec2<f32>(f32(sumX), f32(sumY)) / massF * invPosScale;
  return center;
}

fn shouldUseCell(cellSize: f32, dist: f32, theta: f32) -> bool {
  return cellSize / dist < theta;
}

fn particleRadius(mass: f32) -> f32 {
  return sqrt(abs(mass));
}

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let count = u32(params.particleCount);
  if (gid.x >= count) {
    return;
  }

  var p = particles[gid.x];
  if (p.mass == 0.0) {
    return;
  }
  let levelCount = u32(params.levelCount);

  var force = vec2<f32>(0.0, 0.0);

  for (var level: u32 = 0u; level < levelCount; level = level + 1u) {
    let size = levelInfo[level].x;
    let offset = levelInfo[level].y;
    let cellSize = params.domainSize / f32(size);

    let parentSize = select(1u, levelInfo[level - 1u].x, level > 0u);
    let parentOffset = select(0u, levelInfo[level - 1u].y, level > 0u);

    let totalCells = size * size;
    for (var idx: u32 = 0u; idx < totalCells; idx = idx + 1u) {
      let mass = atomicLoad(&levelMass[offset + idx]);
      if (mass == 0) {
        continue;
      }

      let sumX = atomicLoad(&levelSumX[offset + idx]);
      let sumY = atomicLoad(&levelSumY[offset + idx]);
      let center = cellCenter(mass, sumX, sumY, params.invPosScale);

      // check whether a coarser (parent) cell can be used instead of this finer cell
      let levelCount = u32(params.levelCount);
      let hasParent = (level + 1u) < levelCount;
      if (hasParent) {
        let childX = idx % size;
        let childY = idx / size;
        // parent is the coarser cell at level+1 (fewer cells per side)
        let parentSize = levelInfo[level + 1u].x;
        let parentOffset = levelInfo[level + 1u].y;
        let parentX = childX / 2u;
        let parentY = childY / 2u;
        let parentIdx = parentOffset + parentY * parentSize + parentX;
        let parentMass = atomicLoad(&levelMass[parentIdx]);
        if (parentMass != 0) {
          let parentSumX = atomicLoad(&levelSumX[parentIdx]);
          let parentSumY = atomicLoad(&levelSumY[parentIdx]);
          let parentCenter = cellCenter(parentMass, parentSumX, parentSumY, params.invPosScale);
          let dxParent = wrapDelta(parentCenter.x - p.pos.x);
          let dyParent = wrapDelta(parentCenter.y - p.pos.y);
          var distParent2 = dxParent * dxParent + dyParent * dyParent;
          let minSep = max(params.softening, 1e-6);
          distParent2 = max(distParent2, minSep * minSep);
          let distParent = sqrt(distParent2);
          let parentCellSize = params.domainSize / f32(parentSize);
          if (shouldUseCell(parentCellSize, distParent, params.theta)) {
            // if parent cell should be used, skip processing this finer cell (it will be approximated by parent)
            continue;
          }
        }
      }

        let dx = wrapDelta(center.x - p.pos.x);
        let dy = wrapDelta(center.y - p.pos.y);
        var dist2 = dx * dx + dy * dy;
        let m1 = abs(p.mass);
        let m2raw = f32(mass);
        let m2 = abs(m2raw * params.invMassScale);
        let radius1 = particleRadius(m1);
        let radius2 = particleRadius(m2);
        let mergeT = (radius1 + radius2) * params.mergeRadius;
        if (dist2 < mergeT * mergeT) {
          continue;
        }
        let t1 = select(1.0, -1.0, p.mass < 0.0);
        let t2 = select(1.0, -1.0, m2raw < 0.0);

        let soft = max(params.softening, 1e-4);
        var dist2s = dist2 + soft * soft;
        let invDist = inverseSqrt(dist2s);
        let forceMag = (params.gravity * m1 * m2) * (invDist * invDist);
        let oppSign = t1 * t2 < 0.0;
        let dir = select(1.0, -1.0, oppSign);
        force += vec2<f32>(dx, dy) * (forceMag * dir * invDist);
    }
  }

  let invMassP = 1.0 / max(abs(p.mass), 1e-6);
  let accel = force * invMassP;
  // integrate velocity with a small clamp to avoid runaway speeds when softening=0
  var newVel = (p.vel + accel * params.dt) * params.damping;
  let maxVel = 5.0;
  if (length(newVel) > maxVel) {
    newVel = normalize(newVel) * maxVel;
  }
  p.vel = newVel;
  // open space: don't wrap positions
  p.pos = p.pos + p.vel * params.dt;
  particles[gid.x] = p;
}
`,
///// DIRECT SUM SHADER
  updateParticlesDirect: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read> particlesSrc: array<Particle>;
@group(0) @binding(1) var<storage, read_write> particlesDst: array<Particle>;
@group(0) @binding(2) var<uniform> params: SimParams;

fn particleRadius(mass: f32) -> f32 {
  return sqrt(abs(mass));
}

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let count = u32(params.particleCount);
  if (gid.x >= count) {
    return;
  }

  var p = particlesSrc[gid.x];
  if (p.mass == 0.0) {
    particlesDst[gid.x] = p;
    return;
  }

  var force = vec2<f32>(0.0, 0.0);

  for (var j: u32 = 0u; j < count; j = j + 1u) {
    if (j == gid.x) {
      continue;
    }
    let other = particlesSrc[j];
    if (other.mass == 0.0) {
      continue;
    }
    let t1 = select(1.0, -1.0, p.mass < 0.0);
    let t2 = select(1.0, -1.0, other.mass < 0.0);
    let sameSign = t1 * t2 > 0.0;
    let oppSign = t1 * t2 < 0.0;

    let dx = other.pos.x - p.pos.x;
    let dy = other.pos.y - p.pos.y;
    let soft = max(params.softening, 1e-4);
    var dist2 = dx * dx + dy * dy + soft * soft;
    let mergeT = (particleRadius(p.mass) + particleRadius(other.mass)) * params.mergeRadius;
    if (dist2 < mergeT * mergeT && sameSign) {
      continue;
    }
    let m1 = abs(p.mass);
    let m2 = abs(other.mass);
    let invDist = inverseSqrt(dist2);
    let forceMag = m2 * invDist * invDist;
    let dir = select(1.0, -1.0, oppSign);
    force += vec2<f32>(dx, dy) * (forceMag * dir * invDist);
  }

  let accel = force * params.gravity;
  var newVel = (p.vel + accel * params.dt) * params.damping;
  p.vel = newVel;
  p.pos = p.pos + p.vel * params.dt;
  particlesDst[gid.x] = p;
}
`,
  clearMerge: `
struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read_write> mergeHead: array<atomic<u32>>;
@group(0) @binding(1) var<storage, read_write> mergeNext: array<u32>;
@group(0) @binding(2) var<storage, read_write> mergeTarget: array<u32>;
@group(0) @binding(3) var<storage, read_write> mergeMass: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> mergeSumX: array<atomic<i32>>;
@group(0) @binding(5) var<storage, read_write> mergeSumY: array<atomic<i32>>;
@group(0) @binding(6) var<storage, read_write> mergeSumVx: array<atomic<i32>>;
@group(0) @binding(7) var<storage, read_write> mergeSumVy: array<atomic<i32>>;
@group(0) @binding(8) var<uniform> params: SimParams;

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let cellCount = ${MERGE_GRID}u * ${MERGE_GRID}u;
  let particleCount = u32(params.particleCount);
  if (gid.x < cellCount) {
    atomicStore(&mergeHead[gid.x], 0xffffffffu);
  }
  if (gid.x < particleCount) {
    mergeNext[gid.x] = 0xffffffffu;
    mergeTarget[gid.x] = 0u;
  atomicStore(&mergeMass[gid.x], 0);
  atomicStore(&mergeSumX[gid.x], 0);
  atomicStore(&mergeSumY[gid.x], 0);
  atomicStore(&mergeSumVx[gid.x], 0);
  atomicStore(&mergeSumVy[gid.x], 0);
  }
}
`,
  mergeDeposit: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read_write> mergeHead: array<atomic<u32>>;
@group(0) @binding(1) var<storage, read_write> mergeNext: array<u32>;
@group(0) @binding(2) var<storage, read> particles: array<Particle>;
@group(0) @binding(3) var<uniform> params: SimParams;

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let count = u32(params.particleCount);
  if (gid.x >= count) {
    return;
  }
  let p = particles[gid.x];
  if (p.mass == 0.0) {
    return;
  }

  if (p.pos.x < 0.0 || p.pos.x > 1.0 || p.pos.y < 0.0 || p.pos.y > 1.0) {
    return;
  }

  let size = ${MERGE_GRID}u;
  let fx = p.pos.x * f32(size);
  let fy = p.pos.y * f32(size);
  let x = u32(clamp(floor(fx), 0.0, f32(size - 1)));
  let y = u32(clamp(floor(fy), 0.0, f32(size - 1)));
  let idx = y * size + x;

  let prev = atomicExchange(&mergeHead[idx], gid.x);
  mergeNext[gid.x] = prev;
}
`,
  mergeTargets: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read> particles: array<Particle>;
@group(0) @binding(1) var<storage, read_write> mergeHead: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read> mergeNext: array<u32>;
@group(0) @binding(3) var<storage, read_write> mergeTarget: array<u32>;
@group(0) @binding(4) var<uniform> params: SimParams;

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let count = u32(params.particleCount);
  if (gid.x >= count) {
    return;
  }

  let p = particles[gid.x];
  if (p.mass == 0.0) {
    mergeTarget[gid.x] = 0u;
    return;
  }

  if (p.pos.x < 0.0 || p.pos.x > 1.0 || p.pos.y < 0.0 || p.pos.y > 1.0) {
    mergeTarget[gid.x] = 0u;
    return;
  }

  let size = ${MERGE_GRID}u;
  let fx = p.pos.x * f32(size);
  let fy = p.pos.y * f32(size);
  let x = i32(clamp(floor(fx), 0.0, f32(size - 1)));
  let y = i32(clamp(floor(fy), 0.0, f32(size - 1)));

  let m1 = abs(p.mass);
  let mergeTBase = sqrt(m1) * params.mergeRadius;
  var best = gid.x;
  var bestDist = mergeTBase * mergeTBase;

  for (var dy: i32 = -1; dy <= 1; dy = dy + 1) {
  for (var dx: i32 = -1; dx <= 1; dx = dx + 1) {
  let nx = clamp(x + dx, 0, i32(size) - 1);
  let ny = clamp(y + dy, 0, i32(size) - 1);
  let idx = u32(ny) * size + u32(nx);

      var current = atomicLoad(&mergeHead[idx]);
      var steps: u32 = 0u;
      loop {
        if (current == 0xffffffffu || steps >= 64u) {
          break;
        }

        if (current < gid.x) {
          let other = particles[current];
          if (other.mass != 0.0) {
            // prevent merging with opposite-sign particles
            if (other.mass * p.mass < 0.0) {
              current = mergeNext[current];
              steps = steps + 1u;
              continue;
            }
            var dxp = other.pos.x - p.pos.x;
            var dyp = other.pos.y - p.pos.y;
            // no wrapping for open space
            let dist2 = dxp * dxp + dyp * dyp;
            let m2 = abs(other.mass);
            let mergeT = (sqrt(m1) + sqrt(m2)) * params.mergeRadius;
            if (dist2 < mergeT * mergeT && dist2 < bestDist) {
              bestDist = dist2;
              best = current;
            }
          }
        }

        current = mergeNext[current];
        steps = steps + 1u;
      }
    }
  }

  mergeTarget[gid.x] = best + 1u;
}
`,
  mergeAccumulate: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read> particles: array<Particle>;
@group(0) @binding(1) var<storage, read> mergeTarget: array<u32>;
@group(0) @binding(2) var<storage, read_write> mergeMass: array<atomic<i32>>;
@group(0) @binding(3) var<storage, read_write> mergeSumX: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> mergeSumY: array<atomic<i32>>;
@group(0) @binding(5) var<storage, read_write> mergeSumVx: array<atomic<i32>>;
@group(0) @binding(6) var<storage, read_write> mergeSumVy: array<atomic<i32>>;
@group(0) @binding(7) var<uniform> params: SimParams;

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let count = u32(params.particleCount);
  if (gid.x >= count) {
    return;
  }
  let p = particles[gid.x];
  if (p.mass == 0.0) {
    return;
  }

  let targetIdx = mergeTarget[gid.x];
  if (targetIdx == 0u) {
    return;
  }
  let owner = targetIdx - 1u;

  let massScaled = i32(round(p.mass * params.massScale));
  let sumX = i32(round(p.pos.x * ${MERGE_POS_SCALE.toFixed(2)} * f32(massScaled)));
  let sumY = i32(round(p.pos.y * ${MERGE_POS_SCALE.toFixed(2)} * f32(massScaled)));
  let sumVx = i32(round(p.vel.x * ${MERGE_VEL_SCALE.toFixed(2)} * f32(massScaled)));
  let sumVy = i32(round(p.vel.y * ${MERGE_VEL_SCALE.toFixed(2)} * f32(massScaled)));

  atomicAdd(&mergeMass[owner], massScaled);
  atomicAdd(&mergeSumX[owner], sumX);
  atomicAdd(&mergeSumY[owner], sumY);
  atomicAdd(&mergeSumVx[owner], sumVx);
  atomicAdd(&mergeSumVy[owner], sumVy);
}
`,
  mergeApply: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct SimParams {
  levelCount: f32,
  particleCount: f32,
  domainSize: f32,
  dt: f32,
  theta: f32,
  softening: f32,
  damping: f32,
  gravity: f32,
  massScale: f32,
  invMassScale: f32,
  posScale: f32,
  invPosScale: f32,
  mergeRadius: f32,
  velScale: f32,
  invVelScale: f32,
};

@group(0) @binding(0) var<storage, read_write> particles: array<Particle>;
@group(0) @binding(1) var<storage, read> mergeTarget: array<u32>;
@group(0) @binding(2) var<storage, read_write> mergeMass: array<atomic<i32>>;
@group(0) @binding(3) var<storage, read_write> mergeSumX: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> mergeSumY: array<atomic<i32>>;
@group(0) @binding(5) var<storage, read_write> mergeSumVx: array<atomic<i32>>;
@group(0) @binding(6) var<storage, read_write> mergeSumVy: array<atomic<i32>>;
@group(0) @binding(7) var<uniform> params: SimParams;

fn cellCenter(mass: i32, sumX: i32, sumY: i32, invPosScale: f32) -> vec2<f32> {
  let massF = f32(mass);
  let center = vec2<f32>(f32(sumX), f32(sumY)) / massF * invPosScale;
  return center;
}

@compute @workgroup_size(${WORKGROUP_SIZE_1D})
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let count = u32(params.particleCount);
  if (gid.x >= count) {
    return;
  }

  var p = particles[gid.x];
  if (p.mass == 0.0) {
    return;
  }

  let targetIdx = mergeTarget[gid.x];
  if (targetIdx == 0u) {
    return;
  }
  let owner = targetIdx - 1u;

  let mass = atomicLoad(&mergeMass[owner]);
  if (mass == 0) {
    return;
  }

    if (owner == gid.x) {
      let sumX = atomicLoad(&mergeSumX[owner]);
      let sumY = atomicLoad(&mergeSumY[owner]);
      let sumVx = atomicLoad(&mergeSumVx[owner]);
      let sumVy = atomicLoad(&mergeSumVy[owner]);
    let center = cellCenter(mass, sumX, sumY, ${ (1 / MERGE_POS_SCALE).toFixed(6)});
    // center is already in domain coordinates; write directly (no wrapping)
    p.pos = center;
  p.mass = f32(mass) * params.invMassScale;
  let massF = max(abs(f32(mass)), 1.0);
  let vx = f32(sumVx) / massF * ${ (1 / MERGE_VEL_SCALE).toFixed(6)};
  let vy = f32(sumVy) / massF * ${ (1 / MERGE_VEL_SCALE).toFixed(6)};
    p.vel = vec2<f32>(vx, vy);
  } else {
    p.mass = 0.0;
    p.vel = vec2<f32>(0.0, 0.0);
  }

  particles[gid.x] = p;
}
`,
  render: `
struct Particle {
  pos: vec2<f32>,
  vel: vec2<f32>,
  mass: f32,
  pad: f32,
};

struct ViewParams {
  center: vec2<f32>,
  zoom: f32,
  aspect: f32,
  showPos: f32,
  showNeg: f32,
};

@group(0) @binding(0) var<storage, read> particles: array<Particle>;
@group(0) @binding(1) var<uniform> view: ViewParams;
@group(0) @binding(2) var particleTexture: texture_2d<f32>;
@group(0) @binding(3) var particleSampler: sampler;

struct VSOut {
  @builtin(position) position: vec4<f32>,
  @location(0) mass: f32,
  @location(1) uv: vec2<f32>,
};

const quad: array<vec2<f32>, 6> = array<vec2<f32>, 6>(
  vec2<f32>(-1.0, -1.0),
  vec2<f32>( 1.0, -1.0),
  vec2<f32>(-1.0,  1.0),
  vec2<f32>(-1.0,  1.0),
  vec2<f32>( 1.0, -1.0),
  vec2<f32>( 1.0,  1.0)
);

const quadUv: array<vec2<f32>, 6> = array<vec2<f32>, 6>(
  vec2<f32>(0.0, 1.0),
  vec2<f32>(1.0, 1.0),
  vec2<f32>(0.0, 0.0),
  vec2<f32>(0.0, 0.0),
  vec2<f32>(1.0, 1.0),
  vec2<f32>(1.0, 0.0)
);

@vertex
fn vsMain(@builtin(vertex_index) vid: u32, @builtin(instance_index) iid: u32) -> VSOut {
  let p = particles[iid];
  let radius = sqrt(abs(p.mass));
  let offset = quad[vid] * (0.0002 + radius * 0.001);
  let viewPos = (p.pos - view.center) * view.zoom;
  let clip = vec2<f32>((viewPos.x + offset.x) * 2.0, (viewPos.y + offset.y) * 2.0);
  var out: VSOut;
  out.position = vec4<f32>(clip.x, clip.y * view.aspect, 0.0, 1.0);
  out.mass = p.mass;
  out.uv = quadUv[vid];
  return out;
}

@fragment
fn fsMain(input: VSOut) -> @location(0) vec4<f32> {
  if (input.mass == 0.0) {
    discard;
  }
  // visibility toggles: hide negative/positive particles when their flag is off
  if (input.mass < 0.0 && view.showNeg < 0.5) {
    discard;
  }
  if (input.mass > 0.0 && view.showPos < 0.5) {
    discard;
  }
  let texel = textureSample(particleTexture, particleSampler, input.uv);
  if (texel.a < 0.02) {
    discard;
  }
  let centered = input.uv * 2.0 - vec2<f32>(1.0, 1.0);
  let r2 = dot(centered, centered);
  if (r2 > 1.0) {
    discard;
  }
  let radial = smoothstep(1.0, 0.7, r2);
  let massMag = abs(input.mass);
  let t = clamp((massMag - ${MASS_RANGE[0].toFixed(2)}) / ${(MASS_RANGE[1] - MASS_RANGE[0]).toFixed(2)}, 0.0, 1.0);
  let posColor = mix(vec3<f32>(0.95, 0.25, 0.2), vec3<f32>(1.0, 0.6, 0.4), t);
  let negColor = mix(vec3<f32>(0.2, 0.45, 1.0), vec3<f32>(0.5, 0.8, 1.0), t);
  let color = select(posColor, negColor, input.mass < 0.0);
  return vec4<f32>(color, 0.9 * texel.a * radial);
}
`
};

async function init() {
  if (!navigator.gpu) {
    setStatus('WebGPU not supported in this browser.');
    return;
  }

  const adapter = await navigator.gpu.requestAdapter();
  if (!adapter) {
    setStatus('No compatible WebGPU adapter found.');
    return;
  }

  const device = await adapter.requestDevice();
  const context = canvas.getContext('webgpu');
  const format = navigator.gpu.getPreferredCanvasFormat();

  function resize() {
    const pixelRatio = window.devicePixelRatio || 1;
    canvas.width = Math.floor(canvas.clientWidth * pixelRatio);
    canvas.height = Math.floor(canvas.clientHeight * pixelRatio);
    context.configure({ device, format, alphaMode: 'opaque' });
  }

  resize();
  window.addEventListener('resize', resize);
  setupInteractions();
  setupControls();

  const particleTexture = await loadTexture(
    device,
    'https://assets.babylonjs.com/textures/flare.png'
  );
  const particleSampler = device.createSampler({
    magFilter: 'linear',
    minFilter: 'linear',
    addressModeU: 'clamp-to-edge',
    addressModeV: 'clamp-to-edge'
  });

  const { info: levelInfo, totalCells } = buildLevelInfo();
  const levelInfoBuffer = createBuffer(device, levelInfo, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);

  // particle count is mutable at runtime; create buffers with current PARTICLE_COUNT
  let particleData = createParticles(PARTICLE_COUNT);
  let particleBufferA = createBuffer(device, particleData, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let particleBufferB = createBuffer(device, particleData, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let useBufferA = true;

  const levelMass = createEmptyBuffer(device, totalCells * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  const levelSumX = createEmptyBuffer(device, totalCells * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  const levelSumY = createEmptyBuffer(device, totalCells * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);

  const mergeCells = MERGE_GRID * MERGE_GRID;
  const mergeHead = createEmptyBuffer(device, mergeCells * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let mergeNext = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let mergeTarget = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let mergeMass = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let mergeSumX = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let mergeSumY = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let mergeSumVx = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  let mergeSumVy = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);

  const simParamsBuffer = createBuffer(device, buildSimParams(), GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST);
  const viewParamsBuffer = createEmptyBuffer(device, 32, GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST);

  const pipeline = (code, entryPoint = 'main') =>
    device.createComputePipeline({
      layout: 'auto',
      compute: {
        module: device.createShaderModule({ code }),
        entryPoint
      }
    });

  const clearPipeline = pipeline(shaderSources.clearLevels);
  const depositPipeline = pipeline(shaderSources.deposit);
  const reducePipeline = pipeline(shaderSources.reduceLevels);
  const updatePipeline = pipeline(shaderSources.updateParticles);
  const updateDirectPipeline = pipeline(shaderSources.updateParticlesDirect);
  const clearMergePipeline = pipeline(shaderSources.clearMerge);
  const mergeDepositPipeline = pipeline(shaderSources.mergeDeposit);
  const mergeTargetsPipeline = pipeline(shaderSources.mergeTargets);
  const mergeAccumulatePipeline = pipeline(shaderSources.mergeAccumulate);
  const mergeApplyPipeline = pipeline(shaderSources.mergeApply);

  const renderPipeline = device.createRenderPipeline({
    layout: 'auto',
    vertex: {
      module: device.createShaderModule({ code: shaderSources.render }),
      entryPoint: 'vsMain'
    },
    fragment: {
      module: device.createShaderModule({ code: shaderSources.render }),
      entryPoint: 'fsMain',
      targets: [{
        format,
        blend: {
          color: {
            srcFactor: 'src-alpha',
            dstFactor: 'one-minus-src-alpha',
            operation: 'add'
          },
          alpha: {
            srcFactor: 'one',
            dstFactor: 'one-minus-src-alpha',
            operation: 'add'
          }
        }
      }]
    },
    primitive: {
      topology: 'triangle-list'
    }
  });

  // bind groups are created via a helper so they can be recreated when particle count changes
  let clearBindGroup, reduceBindGroup;
  let particleBindGroupsA, particleBindGroupsB;
  let updateDirectBindGroupAB, updateDirectBindGroupBA;

  function createParticleBindGroups(particleBuffer) {
    const groups = {};

    groups.deposit = device.createBindGroup({
      layout: depositPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: levelMass } },
        { binding: 1, resource: { buffer: levelSumX } },
        { binding: 2, resource: { buffer: levelSumY } },
        { binding: 3, resource: { buffer: particleBuffer } },
        { binding: 4, resource: { buffer: simParamsBuffer } },
        { binding: 5, resource: { buffer: levelInfoBuffer } }
      ]
    });

    groups.update = device.createBindGroup({
      layout: updatePipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: particleBuffer } },
        { binding: 1, resource: { buffer: levelInfoBuffer } },
        { binding: 2, resource: { buffer: levelMass } },
        { binding: 3, resource: { buffer: levelSumX } },
        { binding: 4, resource: { buffer: levelSumY } },
        { binding: 5, resource: { buffer: simParamsBuffer } }
      ]
    });

    groups.render = device.createBindGroup({
      layout: renderPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: particleBuffer } },
        { binding: 1, resource: { buffer: viewParamsBuffer } },
        { binding: 2, resource: particleTexture.createView() },
        { binding: 3, resource: particleSampler }
      ]
    });

    groups.clearMerge = device.createBindGroup({
      layout: clearMergePipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: mergeHead } },
        { binding: 1, resource: { buffer: mergeNext } },
        { binding: 2, resource: { buffer: mergeTarget } },
        { binding: 3, resource: { buffer: mergeMass } },
        { binding: 4, resource: { buffer: mergeSumX } },
        { binding: 5, resource: { buffer: mergeSumY } },
        { binding: 6, resource: { buffer: mergeSumVx } },
        { binding: 7, resource: { buffer: mergeSumVy } },
        { binding: 8, resource: { buffer: simParamsBuffer } }
      ]
    });

    groups.mergeDeposit = device.createBindGroup({
      layout: mergeDepositPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: mergeHead } },
        { binding: 1, resource: { buffer: mergeNext } },
        { binding: 2, resource: { buffer: particleBuffer } },
        { binding: 3, resource: { buffer: simParamsBuffer } }
      ]
    });

    groups.mergeTargets = device.createBindGroup({
      layout: mergeTargetsPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: particleBuffer } },
        { binding: 1, resource: { buffer: mergeHead } },
        { binding: 2, resource: { buffer: mergeNext } },
        { binding: 3, resource: { buffer: mergeTarget } },
        { binding: 4, resource: { buffer: simParamsBuffer } }
      ]
    });

    groups.mergeAccumulate = device.createBindGroup({
      layout: mergeAccumulatePipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: particleBuffer } },
        { binding: 1, resource: { buffer: mergeTarget } },
        { binding: 2, resource: { buffer: mergeMass } },
        { binding: 3, resource: { buffer: mergeSumX } },
        { binding: 4, resource: { buffer: mergeSumY } },
        { binding: 5, resource: { buffer: mergeSumVx } },
        { binding: 6, resource: { buffer: mergeSumVy } },
        { binding: 7, resource: { buffer: simParamsBuffer } }
      ]
    });

    groups.mergeApply = device.createBindGroup({
      layout: mergeApplyPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: particleBuffer } },
        { binding: 1, resource: { buffer: mergeTarget } },
        { binding: 2, resource: { buffer: mergeMass } },
        { binding: 3, resource: { buffer: mergeSumX } },
        { binding: 4, resource: { buffer: mergeSumY } },
        { binding: 5, resource: { buffer: mergeSumVx } },
        { binding: 6, resource: { buffer: mergeSumVy } },
        { binding: 7, resource: { buffer: simParamsBuffer } }
      ]
    });

    return groups;
  }

  function createBindGroups() {
    clearBindGroup = device.createBindGroup({
      layout: clearPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: levelMass } },
        { binding: 1, resource: { buffer: levelSumX } },
        { binding: 2, resource: { buffer: levelSumY } },
        { binding: 3, resource: { buffer: simParamsBuffer } },
        { binding: 4, resource: { buffer: levelInfoBuffer } }
      ]
    });

    reduceBindGroup = device.createBindGroup({
      layout: reducePipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: levelMass } },
        { binding: 1, resource: { buffer: levelSumX } },
        { binding: 2, resource: { buffer: levelSumY } },
        { binding: 3, resource: { buffer: simParamsBuffer } },
        { binding: 4, resource: { buffer: levelInfoBuffer } }
      ]
    });

    particleBindGroupsA = createParticleBindGroups(particleBufferA);
    particleBindGroupsB = createParticleBindGroups(particleBufferB);

    updateDirectBindGroupAB = device.createBindGroup({
      layout: updateDirectPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: particleBufferA } },
        { binding: 1, resource: { buffer: particleBufferB } },
        { binding: 2, resource: { buffer: simParamsBuffer } }
      ]
    });

    updateDirectBindGroupBA = device.createBindGroup({
      layout: updateDirectPipeline.getBindGroupLayout(0),
      entries: [
        { binding: 0, resource: { buffer: particleBufferB } },
        { binding: 1, resource: { buffer: particleBufferA } },
        { binding: 2, resource: { buffer: simParamsBuffer } }
      ]
    });
  }

  // create initial bind groups
  createBindGroups();

  // allow regeneration of particle buffers at runtime
  function regenerateParticles(newCount) {
    const count = Math.max(1, Math.floor(Number(newCount) || PARTICLE_COUNT));
    PARTICLE_COUNT = count;

  // recreate particle buffers with new data
  particleData = createParticles(PARTICLE_COUNT);
  particleBufferA = createBuffer(device, particleData, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  particleBufferB = createBuffer(device, particleData, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
  useBufferA = true;

    // recreate merge-related buffers sized to PARTICLE_COUNT
    mergeNext = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
    mergeTarget = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
    mergeMass = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
    mergeSumX = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
    mergeSumY = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
    mergeSumVx = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);
    mergeSumVy = createEmptyBuffer(device, PARTICLE_COUNT * 4, GPUBufferUsage.STORAGE | GPUBufferUsage.COPY_DST);

  // recreate bind groups to point to new buffers
  createBindGroups();

    setStatus(`Regenerated ${PARTICLE_COUNT} particles`);
  }

  // expose to UI (setupControls registers handler earlier)
  window.regenerateParticles = regenerateParticles;

  // if solver was changed earlier before regenerateParticles existed, perform the deferred reset now
  if (resetOnMethodChange) {
    window.regenerateParticles(PARTICLE_COUNT);
    resetOnMethodChange = false;
    setStatus('Simulation reset after solver change');
  }

  function dispatch1D(pass, pipelineRef, bindGroup, count) {
    pass.setPipeline(pipelineRef);
    pass.setBindGroup(0, bindGroup);
    pass.dispatchWorkgroups(Math.ceil(count / WORKGROUP_SIZE_1D));
  }

  function dispatch2D(pass, pipelineRef, bindGroup, size) {
    pass.setPipeline(pipelineRef);
    pass.setBindGroup(0, bindGroup);
    const groups = Math.ceil(size / WORKGROUP_SIZE);
    pass.dispatchWorkgroups(groups, groups);
  }

  function updateViewParams() {
    const aspect = canvas.height === 0 ? 1 : canvas.width / canvas.height;
    device.queue.writeBuffer(viewParamsBuffer, 0, buildViewParams(aspect));
  }

  updateViewParams();
  window.addEventListener('resize', updateViewParams);


  let frame = 0;
  let lastTime = performance.now();

  function frameLoop() {
    const now = performance.now();
    const dt = (now - lastTime) / 1000;
    lastTime = now;

    updateViewParams();

    const encoder = device.createCommandEncoder();

    // perform compute passes when not paused, or when a single-step was requested
    let currentGroups = useBufferA ? particleBindGroupsA : particleBindGroupsB;
    const shouldCompute = !isPaused || stepRequested;
    if (shouldCompute) {
      const pass = encoder.beginComputePass();

      if (uiState.useDirect) {
        const simParams = buildSimParams();
        device.queue.writeBuffer(simParamsBuffer, 0, simParams);
        const updateGroup = useBufferA ? updateDirectBindGroupAB : updateDirectBindGroupBA;
        dispatch1D(pass, updateDirectPipeline, updateGroup, PARTICLE_COUNT);
        // swap target buffer
        useBufferA = !useBufferA;
        currentGroups = useBufferA ? particleBindGroupsA : particleBindGroupsB;
      } else {
        dispatch1D(pass, clearPipeline, clearBindGroup, totalCells);
        dispatch1D(pass, depositPipeline, currentGroups.deposit, PARTICLE_COUNT);

        for (let level = 1; level < LEVEL_COUNT; level += 1) {
          const levelParams = buildSimParams();
          levelParams[3] = level;
          device.queue.writeBuffer(simParamsBuffer, 0, levelParams);
          dispatch2D(pass, reducePipeline, reduceBindGroup, GRID_LEVELS[level]);
        }

        const simParams = buildSimParams();
        device.queue.writeBuffer(simParamsBuffer, 0, simParams);
        dispatch1D(pass, updatePipeline, currentGroups.update, PARTICLE_COUNT);
      }

      if (uiState.enableMerging) {
        dispatch1D(pass, clearMergePipeline, currentGroups.clearMerge, Math.max(mergeCells, PARTICLE_COUNT));
        dispatch1D(pass, mergeDepositPipeline, currentGroups.mergeDeposit, PARTICLE_COUNT);
        dispatch1D(pass, mergeTargetsPipeline, currentGroups.mergeTargets, PARTICLE_COUNT);
        dispatch1D(pass, mergeAccumulatePipeline, currentGroups.mergeAccumulate, PARTICLE_COUNT);
        dispatch1D(pass, mergeApplyPipeline, currentGroups.mergeApply, PARTICLE_COUNT);
      }
      pass.end();
      // clear single-step request after performing the compute pass
      stepRequested = false;
    }

    const colorAttachment = {
      view: context.getCurrentTexture().createView(),
      // set background to white per user request
      clearValue: { r: 0.0, g: 0.0, b: 0.0, a: 1.0 },
      loadOp: 'clear',
      storeOp: 'store'
    };

    const renderPass = encoder.beginRenderPass({ colorAttachments: [colorAttachment] });
  renderPass.setPipeline(renderPipeline);
  renderPass.setBindGroup(0, currentGroups.render);
    renderPass.draw(6, PARTICLE_COUNT, 0, 0);
    renderPass.end();

    device.queue.submit([encoder.finish()]);

    frame += 1;
    if (frame % 60 === 0) {
      setStatus(`Particles: ${PARTICLE_COUNT}\nLevels: ${LEVEL_COUNT}\nFrame dt: ${(dt * 1000).toFixed(2)} ms`);
    }

    requestAnimationFrame(frameLoop);
  }

  setStatus('Running Barnes-Hut solver on GPU.');
  requestAnimationFrame(frameLoop);
}

init();

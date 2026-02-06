const MAX_PHYSICS_CAPACITY = 100000;
const FLASH_CAPACITY = 120;

const STATE = {
    G: 10.0,
    initV: 0.0,
    dt: 0.6,
    mergeReach: 0.1,
    initialCount: 2000,
    cameraZoom: 0.5,
    showPos: true,
    showNeg: true,
    particleScale: 1.0,
    velocityMode: 'zero'
};

const dataPos = new Float32Array(MAX_PHYSICS_CAPACITY * 2);
const dataVel = new Float32Array(MAX_PHYSICS_CAPACITY * 2);
const dataProps = new Float32Array(MAX_PHYSICS_CAPACITY * 2);
const dataActive = new Uint8Array(MAX_PHYSICS_CAPACITY);
const colorBuffer = new Float32Array(MAX_PHYSICS_CAPACITY * 3);
const sizeBuffer = new Float32Array(MAX_PHYSICS_CAPACITY);

const flashPos = new Float32Array(FLASH_CAPACITY * 3);
const flashLife = new Float32Array(FLASH_CAPACITY);
const flashColor = new Float32Array(FLASH_CAPACITY * 3);
let flashIndex = 0;

let scene, camera, renderer, pointsSystem, geometry, flashSystem, flashGeometry;
let grid = new Map();
let currentActiveCount = 0;
let maxMass = 1.0;

let isDragging = false;
let lastMousePos = { x: 0, y: 0 };
let lastPinchDist = 0;

init();
setupInputs();
animate();

function init() {
    const container = document.getElementById('canvas-container');
    scene = new THREE.Scene();
    const aspect = window.innerWidth / window.innerHeight;
    const frustumSize = 5000;
    camera = new THREE.OrthographicCamera(frustumSize * aspect / -2, frustumSize * aspect / 2, frustumSize / 2, frustumSize / -2, 1, 2000);
    camera.position.z = 100;

    renderer = new THREE.WebGLRenderer({ antialias: false, powerPreference: "high-performance" });
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    container.appendChild(renderer.domElement);

    geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(MAX_PHYSICS_CAPACITY * 3), 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colorBuffer, 3));
    geometry.setAttribute('size', new THREE.BufferAttribute(sizeBuffer, 1));

    const material = new THREE.ShaderMaterial({
        transparent: true, blending: THREE.AdditiveBlending, depthTest: false,
        vertexShader: `
            attribute float size;
            attribute vec3 color;
            varying vec3 vColor;
            void main() {
                vColor = color;
                gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
                gl_PointSize = size;
            }`,
        fragmentShader: `
            varying vec3 vColor;
            void main() {
                float dist = length(gl_PointCoord - vec2(0.5));
                if (dist > 0.5) discard;
                float core = smoothstep(0.5, 0.4, dist) * 0.9 + smoothstep(0.4, 0.0, dist) * 0.1;
                gl_FragColor = vec4(vColor, core);
            }`
    });

    pointsSystem = new THREE.Points(geometry, material);
    scene.add(pointsSystem);

    flashGeometry = new THREE.BufferGeometry();
    flashGeometry.setAttribute('position', new THREE.BufferAttribute(flashPos, 3));
    flashGeometry.setAttribute('color', new THREE.BufferAttribute(flashColor, 3));
    flashGeometry.setAttribute('life', new THREE.BufferAttribute(flashLife, 1));

    const flashMaterial = new THREE.ShaderMaterial({
        transparent: true, blending: THREE.AdditiveBlending, depthTest: false,
        vertexShader: `attribute float life; attribute vec3 color; varying vec3 vColor; varying float vLife; void main() { vColor = color; vLife = life; gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0); gl_PointSize = (30.0 + (1.0 - life) * 100.0); }`,
        fragmentShader: `varying vec3 vColor; varying float vLife; void main() { float dist = length(gl_PointCoord - vec2(0.5)); if (dist > 0.5) discard; float intensity = (1.0 - dist * 2.0) * vLife; gl_FragColor = vec4(vColor * 1.5, intensity); }`
    });

    flashSystem = new THREE.Points(flashGeometry, flashMaterial);
    scene.add(flashSystem);
    resetSim();
}

function updateZoom(factor) {
    STATE.cameraZoom = Math.min(Math.max(STATE.cameraZoom * factor, 0.01), 100);
    camera.zoom = STATE.cameraZoom;
    camera.updateProjectionMatrix();
    document.getElementById('zoom-val').innerText = STATE.cameraZoom.toFixed(2) + 'x';
}

function setupInputs() {
    const inputs = {
        'input-g': ['G', 'g-label', 1],
        'input-v': ['initV', 'iv-label', 1],
        'input-dt': ['dt', 'dt-label', 2],
        'input-merge': ['mergeReach', 'merge-label', 2],
        'input-p': ['initialCount', 'p-label', 0],
        'input-scale': ['particleScale', 'scale-label', 1]
    };

    Object.entries(inputs).forEach(([id, [key, labelId, fixed]]) => {
        const el = document.getElementById(id);
        el.addEventListener('input', (e) => {
            const val = parseFloat(e.target.value);
            STATE[key] = val;
            document.getElementById(labelId).innerText = val.toLocaleString(undefined, { minimumFractionDigits: fixed });
            if (key === 'particleScale') {
                for(let i=0; i<MAX_PHYSICS_CAPACITY; i++) if(dataActive[i]) updateVisual(i);
            }
        });
    });

    document.querySelectorAll('input[name="v-mode"]').forEach(radio => {
        radio.addEventListener('change', (e) => {
            STATE.velocityMode = e.target.value;
        });
    });

    document.getElementById('show-pos').addEventListener('change', (e) => { STATE.showPos = e.target.checked; });
    document.getElementById('show-neg').addEventListener('change', (e) => { STATE.showNeg = e.target.checked; });

    document.getElementById('btn-reset').addEventListener('click', resetSim);

    const uiToggle = document.getElementById('toggle-ui');
    uiToggle.addEventListener('click', () => {
        const panels = document.querySelectorAll('.panel');
        const isHidden = panels[0].classList.toggle('hidden');
        panels[1].classList.toggle('hidden');
        uiToggle.innerText = isHidden ? 'Show UI' : 'Hide UI';
    });

    window.addEventListener('wheel', (e) => { e.preventDefault(); updateZoom(e.deltaY > 0 ? 0.9 : 1.1); }, { passive: false });

    // Drag and Pan interactions
    window.addEventListener('touchstart', (e) => {
        if(e.target.closest('.panel') || e.target.id === 'toggle-ui') return;
        if (e.touches.length === 1) {
            isDragging = true;
            lastMousePos = { x: e.touches[0].clientX, y: e.touches[0].clientY };
        } else if (e.touches.length === 2) {
            lastPinchDist = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
        }
    });

    window.addEventListener('touchmove', (e) => {
        if (e.touches.length === 1 && isDragging) {
            const scale = (1 / STATE.cameraZoom) * (5000 / window.innerHeight);
            camera.position.x -= (e.touches[0].clientX - lastMousePos.x) * scale;
            camera.position.y += (e.touches[0].clientY - lastMousePos.y) * scale;
            lastMousePos = { x: e.touches[0].clientX, y: e.touches[0].clientY };
        } else if (e.touches.length === 2) {
            const dist = Math.hypot(e.touches[0].clientX - e.touches[1].clientX, e.touches[0].clientY - e.touches[1].clientY);
            if (lastPinchDist > 0) updateZoom(dist / lastPinchDist);
            lastPinchDist = dist;
        }
    }, { passive: false });

    window.addEventListener('touchend', () => { isDragging = false; lastPinchDist = 0; });
    window.addEventListener('mousedown', (e) => {
        if(e.target.closest('.panel') || e.target.id === 'toggle-ui') return;
        isDragging = true; lastMousePos = { x: e.clientX, y: e.clientY };
    });
    window.addEventListener('mousemove', (e) => {
        if (!isDragging) return;
        const scale = (1 / STATE.cameraZoom) * (5000 / window.innerHeight);
        camera.position.x -= (e.clientX - lastMousePos.x) * scale;
        camera.position.y += (e.clientY - lastMousePos.y) * scale;
        lastMousePos = { x: e.clientX, y: e.clientY };
    });
    window.addEventListener('mouseup', () => isDragging = false);
    window.addEventListener('resize', onWindowResize);
}

function resetSim() {
    currentActiveCount = STATE.initialCount;
    maxMass = 1.0;
    grid.clear();
    camera.position.set(0, 0, 100);
    STATE.cameraZoom = 0.5;
    camera.zoom = 0.5;
    camera.updateProjectionMatrix();

    const positions = geometry.attributes.position.array;
    for (let i = 0; i < MAX_PHYSICS_CAPACITY; i++) {
        if (i < STATE.initialCount) {
            dataActive[i] = 1;

            const angle = Math.random() * Math.PI * 2;
            const r = Math.sqrt(Math.random()) * 1200;
            const px = Math.cos(angle) * r;
            const py = Math.sin(angle) * r;
            dataPos[i*2] = px;
            dataPos[i*2+1] = py;

            if (STATE.velocityMode === 'orbital') {
                const orbitalSpeed = STATE.initV * (0.8 + Math.random() * 0.4);
                const noise = 0.1;
                dataVel[i*2] = -Math.sin(angle) * orbitalSpeed + (Math.random() - 0.5) * orbitalSpeed * noise;
                dataVel[i*2+1] = Math.cos(angle) * orbitalSpeed + (Math.random() - 0.5) * orbitalSpeed * noise;
            } else {
                dataVel[i*2] = 0;
                dataVel[i*2+1] = 0;
            }

            dataProps[i*2] = 1.0;
            dataProps[i*2+1] = Math.random() > 0.5 ? 1 : -1;
            updateVisual(i);
        } else {
            dataActive[i] = 0;
            positions[i*3] = 999999;
        }
    }
}

function updatePhysics() {
    grid.clear();
    const invGSize = 1.0 / 200;
    for (let i = 0; i < MAX_PHYSICS_CAPACITY; i++) {
        if (dataActive[i] === 0) continue;
        const k = (Math.floor(dataPos[i*2]*invGSize)*31) ^ Math.floor(dataPos[i*2+1]*invGSize);
        if (!grid.has(k)) grid.set(k, []);
        grid.get(k).push(i);
    }

    const posArr = geometry.attributes.position.array;
    for (let i = 0; i < MAX_PHYSICS_CAPACITY; i++) {
        if (dataActive[i] === 0) {
            posArr[i*3] = 999999;
            continue;
        }

        const type = dataProps[i*2+1];
        if ((type > 0 && !STATE.showPos) || (type < 0 && !STATE.showNeg)) {
            posArr[i*3] = 999999;
        } else {
            posArr[i*3] = dataPos[i*2];
            posArr[i*3+1] = dataPos[i*2+1];
        }

        let fx = 0, fy = 0;
        const p1x = dataPos[i*2], p1y = dataPos[i*2+1], m1 = dataProps[i*2], t1 = type;
        const gx = Math.floor(p1x*invGSize), gy = Math.floor(p1y*invGSize);
        const limit = STATE.initialCount > 50000 ? 50 : 200;
        let c = 0;

        outer: for (let dx = -1; dx <= 1; dx++) {
            for (let dy = -1; dy <= 1; dy++) {
                const cell = grid.get(((gx+dx)*31)^(gy+dy));
                if (!cell) continue;
                for (let j of cell) {
                    if (i === j || dataActive[j] === 0) continue;
                    const dxP = dataPos[j*2]-p1x, dyP = dataPos[j*2+1]-p1y, d2 = dxP*dxP+dyP*dyP;
                    const mergeT = (Math.sqrt(m1)+Math.sqrt(dataProps[j*2]))*1.8 * STATE.mergeReach;

                    if (d2 < mergeT*mergeT) {
                        if (m1 >= dataProps[j*2]) {
                            const nm = m1+dataProps[j*2];
                            // ONLY FLASH FOR POSITIVE PARTICLES (Red)
                            if(type > 0 && STATE.initialCount < 40000) {
                                // WHITE FLASH
                                createFlash(dataPos[j*2], dataPos[j*2+1], 1.0, 1.0, 1.0);
                            }
                            dataVel[i*2] = (dataVel[i*2]*m1 + dataVel[j*2]*dataProps[j*2])/nm;
                            dataVel[i*2+1] = (dataVel[i*2+1]*m1 + dataVel[j*2+1]*dataProps[j*2])/nm;
                            dataPos[i*2] = (p1x*m1+dataPos[j*2]*dataProps[j*2])/nm;
                            dataPos[i*2+1] = (p1y*m1+dataPos[j*2+1]*dataProps[j*2])/nm;
                            dataProps[i*2] = nm; dataActive[j] = 0; currentActiveCount--;
                            if (nm > maxMass) maxMass = nm;
                            updateVisual(i);
                        }
                        continue;
                    }
                    const d = Math.sqrt(d2+100), force = (STATE.G*m1*dataProps[j*2])/(d2+400);
                    const dir = (t1*dataProps[j*2+1]) > 0 ? 1 : -2.5;
                    fx += force*(dxP/d)*dir; fy += force*(dyP/d)*dir;
                    if (++c > limit) break outer;
                }
            }
        }
        dataVel[i*2] += (fx/m1)*STATE.dt; dataVel[i*2+1] += (fy/m1)*STATE.dt;
    }

    for (let i = 0; i < MAX_PHYSICS_CAPACITY; i++) {
        if (dataActive[i] === 1) {
            dataPos[i*2] += dataVel[i*2]*STATE.dt; dataPos[i*2+1] += dataVel[i*2+1]*STATE.dt;
        }
    }
    geometry.attributes.position.needsUpdate = true;
}

function updateVisual(i) {
    const m = dataProps[i*2], type = dataProps[i*2+1];
    const baseSize = (STATE.initialCount > 30000 ? 1 : 2);
    sizeBuffer[i] = (baseSize + Math.pow(m, 0.45)*5) * STATE.particleScale;

    const colors = geometry.attributes.color.array;
    if (type > 0) {
        colors[i*3] = 1.0; colors[i*3+1] = 0.2; colors[i*3+2] = 0.1;
    } else {
        colors[i*3] = 0.1; colors[i*3+1] = 0.5; colors[i*3+2] = 1.0;
    }
    geometry.attributes.color.needsUpdate = true; geometry.attributes.size.needsUpdate = true;
}

function updateFlashes() {
    for (let i = 0; i < FLASH_CAPACITY; i++) {
        if (flashLife[i] > 0) {
            flashLife[i] -= 0.05;
            if (flashLife[i] < 0) flashLife[i] = 0;
        } else {
            flashPos[i*3] = 999999;
        }
    }
    flashGeometry.attributes.position.needsUpdate = true;
    flashGeometry.attributes.life.needsUpdate = true;
}

function animate() {
    requestAnimationFrame(animate);
    updatePhysics(); updateFlashes();
    document.getElementById('count').innerText = currentActiveCount.toLocaleString();
    document.getElementById('max-mass').innerText = maxMass.toFixed(1);
    renderer.render(scene, camera);
}

function createFlash(x, y, r, g, b) {
    const idx = flashIndex % FLASH_CAPACITY;
    flashPos[idx*3] = x; flashPos[idx*3+1] = y; flashPos[idx*3+2] = 2;
    flashColor[idx*3] = r; flashColor[idx*3+1] = g; flashColor[idx*3+2] = b;
    flashLife[idx] = 1.0; flashIndex++;
}

function onWindowResize() {
    const aspect = window.innerWidth / window.innerHeight;
    const frustumSize = 5000;
    camera.left = -frustumSize * aspect / 2; camera.right = frustumSize * aspect / 2;
    camera.top = frustumSize / 2; camera.bottom = -frustumSize / 2;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

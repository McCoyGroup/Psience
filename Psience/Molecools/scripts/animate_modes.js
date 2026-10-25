// Controls for the actual sampled scans produced by animate_mode/animate_coordinate.
const data = component.data;
const nodes = component.nodes;
const controls = Object.fromEntries(["mode", "scale", "label"].map(
    role => [role, component.query('[data-role="' + role + '"]')]
));
const byId = new Map(Array.from(component.scene.querySelectorAll('[id]'), node => [node.id, node]));
const state = {mode: 0, scale: 1, fraction: 0, running: data.running};

const scans = data.modes.map(mode => displacement_scale.capture(mode, data.allowAdjustments, byId));

function segment(node, start, end, radius, radialScale = [1, 1]) {
    const delta = end.map((value, axis) => value - start[axis]);
    const length = Math.hypot(...delta);
    const horizontal = Math.hypot(delta[0], delta[2]);
    let rotation;
    if (length < 1e-12 || (horizontal < 1e-12 && delta[1] >= 0)) {
        rotation = [0, 0, 1, 0];
    } else if (horizontal < 1e-12) {
        rotation = [1, 0, 0, Math.PI];
    } else {
        rotation = [delta[2] / horizontal, 0, -delta[0] / horizontal,
                    Math.acos(Math.max(-1, Math.min(1, delta[1] / length)))];
    }
    component.setField(node, "translation", start.map((value, axis) => (value + end[axis]) / 2));
    component.setField(node, "rotation", rotation);
    component.setField(node, "scale", length < 1e-12 ? [0, 0, 0] : [radius * radialScale[0], length, radius * radialScale[1]]);
}

function drawBindings(bindings, points) {
    for (const binding of bindings) {
        if (binding.kind === "position") {
            component.setField(binding.node, "translation", points[binding.point]);
        } else if (binding.kind === "segment") {
            segment(binding.node, points[binding.start], points[binding.end], binding.radius, binding.radialScale);
        } else if (binding.kind === "coordinates") {
            component.setField(binding.node, "point", points.slice(binding.point, binding.point + binding.count));
        }
    }
}

function draw(scan, points) {
    // Playback changes only the molecule's bindings. Motion arrows are static
    // and are updated independently by displacement_scale on a scale/mode edit.
    drawBindings(scan.bindings, points);
}

function pointsFromEvent(value, count) {
    if (typeof value === "string") value = value.trim().split(/[\s,]+/).map(Number);
    if (!value || !value.length) return null;
    let points;
    if (typeof value[0] === "number") {
        if (value.length !== count * 3) return null;
        points = Array.from({length: count}, (_, i) => Array.from(value.slice(3 * i, 3 * i + 3)));
    } else {
        points = Array.from(value, point => "x" in point ? [point.x, point.y, point.z] : Array.from(point));
    }
    return points.length === count && points.every(point => point.length === 3 && point.every(Number.isFinite)) ? points : null;
}

function drawFraction(scan, fraction) {
    let upper = scan.keys.findIndex(key => key >= fraction);
    if (upper < 0) upper = scan.keys.length - 1;
    const lower = Math.max(0, upper - 1);
    const interval = scan.keys[upper] - scan.keys[lower];
    const weight = interval > 0 ? (fraction - scan.keys[lower]) / interval : 0;
    draw(scan, scan.frames[lower].map((point, atom) => point.map((value, axis) =>
        value + weight * (scan.frames[upper][atom][axis] - value)
    )));
}

function getState() {
    return {...state, ...playback.getState(), modeIndex: scans[state.mode].index};
}

function syncControls() {
    if (controls.mode) controls.mode.value = String(state.mode);
    if (controls.scale) controls.scale.value = String(state.scale);
    if (controls.label) controls.label.textContent = scans[state.mode].label;
    component.root.dispatchEvent(new CustomEvent("scanchange", {detail: getState(), bubbles: true}));
}

function rebuild() {
    const scan = scans[state.mode];
    const interpolator = nodes["interpolator" + state.mode];
    scan.frames = displacement_scale.frames(scan, state.scale);
    component.setField(interpolator, "key", scan.keys);
    component.setField(interpolator, "keyValue", scan.frames);
    displacement_scale.updateArrows(scan, state.scale);
    drawFraction(scan, state.fraction);
    component.setField(interpolator, "set_fraction", state.fraction);
    syncControls();
}

function setMode(index) {
    if (!Number.isInteger(index) || index < 0 || index >= scans.length) throw new RangeError("Mode index is outside this animation");
    state.mode = index;
    component.setField(nodes.switch, "whichChoice", index);
    rebuild();
}

function selectMode(sourceIndex) {
    const index = scans.findIndex(scan => scan.index === sourceIndex);
    if (index < 0) throw new RangeError("Mode is not included in this animation");
    setMode(index);
}

function setScale(scale) {
    if (!data.allowAdjustments) throw new Error("This animation does not allow adjustments");
    if (!Number.isFinite(scale) || scale < 0) throw new RangeError("Displacement scale must be finite and nonnegative");
    state.scale = scale;
    rebuild();
}

scans.forEach((scan, index) => component.on(nodes["interpolator" + index], "outputchange", event => {
    if (index !== state.mode || event.fieldName !== "value_changed") return;
    const points = pointsFromEvent(event.value, scan.reference.length);
    if (points) draw(scan, points);
}));
if (controls.mode) component.on(controls.mode, "change", () => setMode(Number(controls.mode.value)));
if (controls.scale) component.on(controls.scale, "input", () => setScale(Number(controls.scale.value)));
const playback = X3DPlayback(component, {
    clock: nodes.clock,
    targets: scans.map((_, index) => nodes["interpolator" + index]),
    duration: data.duration, running: data.running, loop: true
});
playback.onChange(snapshot => {
    state.fraction = snapshot.fraction;
    state.running = snapshot.running;
});
component.setField(nodes.switch, "whichChoice", 0);
rebuild();
return {...playback, setMode, selectMode, setScale, getState};

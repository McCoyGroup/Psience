// This script owns displacement scaling, including the independent, static
// motion-envelope arrows. Arrow points never belong to a time interpolator.
const displacement_scale = (() => {
    const copyPoints = points => points.map(point => point.slice());
    function bindings(values, byId) {
        return values.map(binding => {
            const node = byId.get(binding.id);
            if (!node) throw new Error("Missing scan animation node: " + binding.id);
            return {...binding, node};
        });
    }
    function capture(mode, allowAdjustments, byId) {
        const original = mode.frames.map(copyPoints);
        const reference = copyPoints(mode.reference);
        const differences = allowAdjustments ? original.map((frame, index) =>
            frame.map((point, atom) => point.map((value, axis) =>
                value - (index === 0 ? reference[atom][axis] : original[index - 1][atom][axis])
            ))
        ) : null;
        const arrowReference = copyPoints(mode.arrowReference || []);
        const arrowAnchor = copyPoints(mode.arrowAnchor || []);
        if (arrowReference.length !== arrowAnchor.length) throw new Error("Motion-arrow anchors do not match their points");
        return {...mode, original, reference, differences, frames: original,
            bindings: bindings(mode.bindings, byId),
            arrowReference, arrowAnchor, arrowBindings: bindings(mode.arrowBindings || [], byId)};
    }
    function frames(scan, scale) {
        if (!scan.differences) return scan.original;
        let current = copyPoints(scan.reference);
        return scan.differences.map(increment => {
            current = current.map((point, atom) => point.map((value, axis) =>
                value + scale * increment[atom][axis]
            ));
            return current;
        });
    }
    function updateArrows(scan, scale) {
        const points = scan.arrowReference.map((point, index) => point.map((value, axis) =>
            scan.arrowAnchor[index][axis] + scale * (value - scan.arrowAnchor[index][axis])
        ));
        drawBindings(scan.arrowBindings, points);
    }
    return {capture, frames, updateArrows};
})();

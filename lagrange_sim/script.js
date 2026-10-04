// ============================================================
// CONSTANTS
// ============================================================

const G = 6.67e-11;


// ============================================================
// INITIAL VALUES
// ============================================================

const M1_initial = 6e24;
const M2_initial = 7e22;
const d = 3.84e8;


// ============================================================
// CANVAS
// ============================================================

const canvas = document.getElementById("simulation");
const ctx = canvas.getContext("2d");


// ============================================================
// SLIDERS
// ============================================================

const M1_slider = document.getElementById("mass1");
const M2_slider = document.getElementById("mass2");

const M1_value = document.getElementById("mass1Value");
const M2_value = document.getElementById("mass2Value");


// ============================================================
// RESIZE
// ============================================================

function resizeCanvas() {

    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;

    update();
}

window.addEventListener("resize", resizeCanvas);


// ============================================================
// VECTOR DISTANCE
// ============================================================

function distance(x1, y1, x2, y2) {

    return Math.sqrt(
        (x1 - x2) ** 2 +
        (y1 - y2) ** 2
    );
}


// ============================================================
// NEWTON'S METHOD FOR L1, L2, L3
// ============================================================

function findCollinearPoint(
    x0,
    x1,
    x2,
    M1,
    M2,
    omega
) {

    function f(x) {

        const r1 = Math.abs(x - x1);
        const r2 = Math.abs(x - x2);

        return (
            G * M1
            * (x - x1)
            / r1 ** 3

            +

            G * M2
            * (x - x2)
            / r2 ** 3

            -

            omega ** 2 * x
        );
    }


    function df(x) {

        const r1 = Math.abs(x - x1);
        const r2 = Math.abs(x - x2);

        return (
            -2 * G * M1 / r1 ** 3

            -

            2 * G * M2 / r2 ** 3

            -

            omega ** 2
        );
    }


    let x = x0;


    for (let i = 0; i < 100; i++) {

        const fValue = f(x);
        const dfValue = df(x);

        const xNew =
            x - fValue / dfValue;


        if (Math.abs(xNew - x) < 1e-3) {

            return xNew;
        }

        x = xNew;
    }


    return x;
}


// ============================================================
// 2-D NEWTON METHOD
// ============================================================

function findLagrangeTriangularPoints(
    x1,
    x2,
    M1,
    M2,
    omega,
    x0,
    y0
) {

    function F(x, y) {

        const r1 = Math.sqrt(
            (x - x1) ** 2 +
            y ** 2
        );

        const r2 = Math.sqrt(
            (x - x2) ** 2 +
            y ** 2
        );


        const dVdx =

            G * M1
            * (x - x1)
            / r1 ** 3

            +

            G * M2
            * (x - x2)
            / r2 ** 3

            -

            omega ** 2 * x;


        const dVdy =

            G * M1
            * y
            / r1 ** 3

            +

            G * M2
            * y
            / r2 ** 3

            -

            omega ** 2 * y;


        return [dVdx, dVdy];
    }


    function J(x, y) {

        const r1 = Math.sqrt(
            (x - x1) ** 2 +
            y ** 2
        );

        const r2 = Math.sqrt(
            (x - x2) ** 2 +
            y ** 2
        );


        const Vxx =

            G * M1 *
            (
                1 / r1 ** 3
                -
                3 * (x - x1) ** 2 / r1 ** 5
            )

            +

            G * M2 *
            (
                1 / r2 ** 3
                -
                3 * (x - x2) ** 2 / r2 ** 5
            )

            -

            omega ** 2;


        const Vxy =

            -3 * G * M1
            * (x - x1)
            * y
            / r1 ** 5

            -

            3 * G * M2
            * (x - x2)
            * y
            / r2 ** 5;


        const Vyy =

            G * M1 *
            (
                1 / r1 ** 3
                -
                3 * y ** 2 / r1 ** 5
            )

            +

            G * M2 *
            (
                1 / r2 ** 3
                -
                3 * y ** 2 / r2 ** 5
            )

            -

            omega ** 2;


        return [
            [Vxx, Vxy],
            [Vxy, Vyy]
        ];
    }


    let x = x0;
    let y = y0;


    for (let i = 0; i < 100; i++) {

        const f = F(x, y);

        const jac = J(x, y);


        const a = jac[0][0];
        const b = jac[0][1];

        const c = jac[1][0];
        const d = jac[1][1];


        const determinant =
            a * d - b * c;


        if (Math.abs(determinant) < 1e-30) {

            break;
        }


        // Solve:
        //
        // J * delta = F
        //

        const dx =
            (
                d * f[0]
                -
                b * f[1]
            )
            / determinant;


        const dy =
            (
                -c * f[0]
                +
                a * f[1]
            )
            / determinant;


        const newX = x - dx;
        const newY = y - dy;


        const difference = Math.sqrt(
            (newX - x) ** 2 +
            (newY - y) ** 2
        );


        x = newX;
        y = newY;


        if (difference < 1e-6) {

            break;
        }
    }


    return [x, y];
}


// ============================================================
// POTENTIAL
// ============================================================

function potential(
    x,
    y,
    x1,
    x2,
    M1,
    M2,
    omega
) {

    const r1 = Math.sqrt(
        (x - x1) ** 2 +
        y ** 2
    );

    const r2 = Math.sqrt(
        (x - x2) ** 2 +
        y ** 2
    );


    if (r1 < 1e7 || r2 < 1e7) {

        return NaN;
    }


    return (

        -G * M1 / r1

        -

        G * M2 / r2

        -

        0.5
        * omega ** 2
        * (x ** 2 + y ** 2)

    );
}


// ============================================================
// COLOR FOR POTENTIAL
// ============================================================

function potentialColor(
    value,
    min,
    max
) {

    let t =
        (value - min)
        /
        (max - min);


    t = Math.max(
        0,
        Math.min(1, t)
    );


    // Inferno-like gradient

    const r = Math.floor(
        255 * Math.min(1, t * 2)
    );

    const g = Math.floor(
        150 * Math.max(0, t - 0.3)
    );

    const b = Math.floor(
        100 * Math.max(0, 0.4 - t)
    );


    return `rgb(${r}, ${g}, ${b})`;
}


// ============================================================
// DRAW POTENTIAL CONTOURS
// ============================================================

function drawPotential(
    x1,
    x2,
    M1,
    M2,
    omega,
    scale
) {

    const gridSize = 100;

    const width = canvas.width;
    const height = canvas.height;

    const L = 1.5 * d;


    const values = [];

    let minV = Infinity;
    let maxV = -Infinity;


    // --------------------------------------------------------
    // Calculate potential grid
    // --------------------------------------------------------

    for (let j = 0; j < gridSize; j++) {

        values[j] = [];


        const y =
            L -
            (2 * L * j)
            / (gridSize - 1);


        for (let i = 0; i < gridSize; i++) {

            const x =
                -L +
                (2 * L * i)
                / (gridSize - 1);


            const V = potential(
                x,
                y,
                x1,
                x2,
                M1,
                M2,
                omega
            );


            values[j][i] = V;


            if (!Number.isNaN(V)) {

                minV =
                    Math.min(minV, V);

                maxV =
                    Math.max(maxV, V);
            }
        }
    }


    // --------------------------------------------------------
    // Draw potential field
    // --------------------------------------------------------

    for (let j = 0; j < gridSize - 1; j++) {

        for (let i = 0; i < gridSize - 1; i++) {

            const V = values[j][i];


            if (Number.isNaN(V)) {
                continue;
            }


            const x =
                -L +
                (2 * L * i)
                / (gridSize - 1);


            const y =
                L -
                (2 * L * j)
                / (gridSize - 1);


            const nextX =
                -L +
                (2 * L * (i + 1))
                / (gridSize - 1);


            const nextY =
                L -
                (2 * L * (j + 1))
                / (gridSize - 1);


            const px =
                canvas.width / 2
                +
                x * scale;


            const py =
                canvas.height / 2
                -
                y * scale;


            const cellWidth =
                (nextX - x) * scale;


            const cellHeight =
                (y - nextY) * scale;


            ctx.fillStyle =
                potentialColor(
                    V,
                    minV,
                    maxV
                );


            ctx.globalAlpha = 0.18;


            ctx.fillRect(
                px,
                py - cellHeight,
                cellWidth,
                cellHeight
            );
        }
    }


    ctx.globalAlpha = 1;
}


// ============================================================
// WORLD → SCREEN
// ============================================================

function worldToScreen(
    x,
    y,
    scale
) {

    return {

        x:
            canvas.width / 2
            +
            x * scale,

        y:
            canvas.height / 2
            -
            y * scale
    };
}


// ============================================================
// DRAW BODY
// ============================================================

function drawBody(
    x,
    y,
    radius,
    color,
    scale
) {

    const p =
        worldToScreen(
            x,
            y,
            scale
        );


    ctx.beginPath();

    ctx.arc(
        p.x,
        p.y,
        Math.max(
            radius * scale,
            5
        ),
        0,
        2 * Math.PI
    );

    ctx.fillStyle = color;

    ctx.fill();
}


// ============================================================
// DRAW LAGRANGE POINT
// ============================================================

function drawLagrangePoint(
    x,
    y,
    label,
    scale
) {

    const p =
        worldToScreen(
            x,
            y,
            scale
        );


    ctx.beginPath();

    ctx.arc(
        p.x,
        p.y,
        4,
        0,
        2 * Math.PI
    );


    ctx.fillStyle = "white";

    ctx.fill();


    ctx.fillStyle = "#11FA00";

    ctx.font = "16px Arial";

    ctx.fillText(
        label,
        p.x + 8,
        p.y + 5
    );
}


// ============================================================
// UPDATE SIMULATION
// ============================================================

function update() {

    const M1 =
        Number(M1_slider.value);

    const M2 =
        Number(M2_slider.value);


    // --------------------------------------------------------
    // Angular velocity
    // --------------------------------------------------------

    const omega =
        Math.sqrt(
            G * (M1 + M2)
            / d ** 3
        );


    // --------------------------------------------------------
    // COM-centered body positions
    // --------------------------------------------------------

    const x1 =
        -M2
        /
        (M1 + M2)
        *
        d;


    const x2 =
        M1
        /
        (M1 + M2)
        *
        d;


    // --------------------------------------------------------
    // Update labels
    // --------------------------------------------------------

    M1_value.textContent =
        M1.toExponential(2);

    M2_value.textContent =
        M2.toExponential(2);


    // --------------------------------------------------------
    // Clear canvas
    // --------------------------------------------------------

    ctx.fillStyle = "black";

    ctx.fillRect(
        0,
        0,
        canvas.width,
        canvas.height
    );


    // --------------------------------------------------------
    // Scale
    // --------------------------------------------------------

    const L = 1.5 * d;

    const scale =
        Math.min(
            canvas.width,
            canvas.height
        )
        /
        (2 * L);


    // --------------------------------------------------------
    // Potential
    // --------------------------------------------------------

    drawPotential(
        x1,
        x2,
        M1,
        M2,
        omega,
        scale
    );


    // ========================================================
    // L1
    // ========================================================

    const L1 =
        findCollinearPoint(
            x1 + 0.85 * d,
            x1,
            x2,
            M1,
            M2,
            omega
        );


    // ========================================================
    // L2
    // ========================================================

    const L2 =
        findCollinearPoint(
            x2 + 0.15 * d,
            x1,
            x2,
            M1,
            M2,
            omega
        );


    // ========================================================
    // L3
    // ========================================================

    const L3 =
        findCollinearPoint(
            x1 - 0.5 * d,
            x1,
            x2,
            M1,
            M2,
            omega
        );


    // ========================================================
    // L4 / L5
    // ========================================================

    const xGuess =
        0.5 * (x1 + x2);

    const yGuess =
        Math.sqrt(3)
        / 2
        * d;


    const L4 =
        findLagrangeTriangularPoints(
            x1,
            x2,
            M1,
            M2,
            omega,
            xGuess,
            yGuess
        );


    const L5 =
        findLagrangeTriangularPoints(
            x1,
            x2,
            M1,
            M2,
            omega,
            xGuess,
            -yGuess
        );


    // ========================================================
    // DRAW BODIES
    // ========================================================

    drawBody(
        x1,
        0,
        6.3e6,
        "#0627B8",
        scale
    );


    drawBody(
        x2,
        0,
        6.3e6 / 6,
        "#B0B0B0",
        scale
    );


    // ========================================================
    // DRAW CENTER OF MASS
    // ========================================================

    const com =
        worldToScreen(
            0,
            0,
            scale
        );


    ctx.beginPath();

    ctx.arc(
        com.x,
        com.y,
        4,
        0,
        2 * Math.PI
    );

    ctx.fillStyle = "red";

    ctx.fill();


    // ========================================================
    // DRAW LAGRANGE POINTS
    // ========================================================

    drawLagrangePoint(
        L1,
        0,
        "L1",
        scale
    );


    drawLagrangePoint(
        L2,
        0,
        "L2",
        scale
    );


    drawLagrangePoint(
        L3,
        0,
        "L3",
        scale
    );


    drawLagrangePoint(
        L4[0],
        L4[1],
        "L4",
        scale
    );


    drawLagrangePoint(
        L5[0],
        L5[1],
        "L5",
        scale
    );


    // ========================================================
    // TITLE
    // ========================================================

    ctx.fillStyle = "white";

    ctx.font = "20px Arial";

    ctx.textAlign = "center";

    ctx.fillText(
        "A potential map of the 2-body system",
        canvas.width / 2,
        35
    );


    ctx.textAlign = "left";
}


// ============================================================
// SLIDER EVENTS
// ============================================================

M1_slider.addEventListener(
    "input",
    update
);

M2_slider.addEventListener(
    "input",
    update
);


// ============================================================
// START
// ============================================================

resizeCanvas();
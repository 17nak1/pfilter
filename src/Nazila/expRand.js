
expRand = {}

expRand.exp_rand = function (unif_rand) {
    var q = [
        0.6931471805599453,
        0.9333736875190459,
        0.9888777961838675,
        0.9984959252914960,
        0.9998292811061389,
        0.9999833164100727,
        0.9999985691438767,
        0.9999998906925558,
        0.9999999924734159,
        0.9999999995283275,
        0.9999999999728814,
        0.9999999999985598,
        0.9999999999999289,
        0.9999999999999968,
        0.9999999999999999,
        1.0000000000000000
    ];
    var a = 0.;
    var u = unif_rand();
    while (u <= 0. || u >= 1.)
        u = unif_rand();
    while (true) {
        u += u;
        if (u > 1.)
            break;
        a += q[0];
    }
    u -= 1.;
    if (u <= q[0])
        return a + u;
    var i = 0;
    var ustar = unif_rand();
    var umin = ustar;
    do {
        ustar = unif_rand();
        if (umin > ustar)
            umin = ustar;
        i++;
    } while (u > q[i]);
    return a + umin * q[0];
}

module.exports = expRand
console.log(Math.random)
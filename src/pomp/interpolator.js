interpolator = function(points) {
    var first, interpolated, leftExtrapolated, rightExtrapolated
    let n = points.length - 1

    if (points.length === 0) {
        return function() {
            return 0
        }
    }
    if (points.length === 1) {
        return function() {
            return points[0][1]
        }
    }

    // points = points.sort(function(a, b) {
    //     return a[0] - b[0]
    // })
    first = points[0]

    interpolated = function(x, a, b) {
        return { pop : a['pop'] + (x - a['time']) * (b['pop'] - a['pop']) / (b['time'] - a['time']),
            birthrate : a['birthrate'] + (x - a['time']) * (b['birthrate'] - a['birthrate']) / (b['time'] - a['time'])
        }
    }

    leftExtrapolated = function(x) {
        var a = points[0]
        var b = points[1]

        return { pop : a['pop'] + (x - a['time']) * (b['pop'] - a['pop']) / (b['time'] - a['time']),
            birthrate : a['birthrate'] + (x - a['time']) * (b['birthrate'] - a['birthrate']) / (b['time'] - a['time'])
        }
    }

    rightExtrapolated = function(x) {
        var a = points[n - 1]
        var b = points[n]

        return { pop : b['pop'] + (x - b['time']) * (b['pop'] - a['pop']) / (b['time'] - a['time']),
            birthrate : b['birthrate'] + (x - b['time']) * (b['birthrate'] - a['birthrate']) / (b['time'] - a['time'])
        }
    }

    return function(x) {
        if (x <= first['time']) {
            return leftExtrapolated(x)
        }
        for (var i = 0; i < n; i += 1) {
            if (x > points[i]['time'] && x <= points[i + 1]['time']) {
                return interpolated(x, points[i], points[i + 1])
            }
        }
        return rightExtrapolated(x)
    }
}
    module.exports = interpolator;
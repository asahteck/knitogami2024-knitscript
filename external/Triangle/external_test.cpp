#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <math.h>
#include <iostream>
#include <vector>
#include <cassert>

/* The vertex data structure.  Each vertex is actually an array of REALs.    */
/*   The number of REALs is unknown until runtime.  An integer boundary      */
/*   marker, and sometimes a pointer to a triangle, is appended after the    */
/*   REALs.                                                                  */
typedef REAL *vertex;


REAL distanceToSegment(REAL x, REAL y, REAL xi, REAL yi, REAL xf, REAL yf) {
    // https://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    REAL a = x - xi;
    REAL b = y - yi;
    REAL c = xf - xi;
    REAL d = yf - yi;

    REAL dot = a * c + b * d;
    REAL len_sq = c * c + d * d;
    REAL temp = -1.0;
    if (len_sq != 0.0) { temp = dot / len_sq; }

    // Find the location on the segment that the point is closest to.
    REAL xx;
    REAL yy;

    if (temp < 0.0) {
        xx = xi;
        yy = yi;
    }
    else if (temp > 1.0) {
        xx = xf;
        yy = yf;
    }
    else {
        xx = xi + temp * c;
        yy = yi + temp * d;
    }

    REAL dx = x - xx;
    REAL dy = y - yy;
    REAL dist = sqrt(dx * dx + dy * dy);

    return dist;
}


struct PointTest {
    REAL x0;
    REAL y0;
    REAL maxRes;
    REAL minRes;
    REAL decayRes;
    REAL offset;

    PointTest(REAL maxRes, REAL minRes, REAL decayRes, REAL x0, REAL y0) :
    maxRes(maxRes),
    minRes(minRes),
    decayRes(decayRes),
    x0(x0), y0(y0),
    offset(0.0)
    {}

    PointTest(REAL maxRes, REAL minRes, REAL decayRes, REAL x0, REAL y0, REAL offset) :
    maxRes(maxRes),
    minRes(minRes),
    decayRes(decayRes),
    x0(x0), y0(y0),
    offset(offset)
    {}

    REAL operator()(REAL x, REAL y) {
        REAL dx = x - x0;
        REAL dy = y - y0;
        REAL dist = sqrt(dx * dx + dy * dy);

        if (dist > offset) {
            dist = dist - offset;
        }
        else {
            return minRes;
        }

        return maxRes - (maxRes - minRes) * exp(-dist / decayRes);
    }
};


struct LineTest {
    REAL maxRes;
    REAL minRes;
    REAL decayRes;
    REAL xi;
    REAL yi;
    REAL xf;
    REAL yf;
    REAL offset;

    LineTest(REAL maxRes, REAL minRes, REAL decayRes, REAL xi, REAL yi, REAL xf, REAL yf) :
    maxRes(maxRes),
    minRes(minRes),
    decayRes(decayRes),
    xi(xi), yi(yi), xf(xf), yf(yf),
    offset(0.0)
    {}

    LineTest(REAL maxRes, REAL minRes, REAL decayRes, REAL xi, REAL yi, REAL xf, REAL yf, REAL offset) :
    maxRes(maxRes),
    minRes(minRes),
    decayRes(decayRes),
    xi(xi), yi(yi), xf(xf), yf(yf),
    offset(offset)
    {}

    REAL operator()(REAL x, REAL y) {
        REAL dist = distanceToSegment(x, y, xi, yi, xf, yf);

        if (dist > offset) {
            dist = dist - offset;
        }
        else {
            return minRes;
        }

        return maxRes - (maxRes - minRes) * exp(-dist / decayRes);
    }
};


struct PolygonTest {
    int nVertices;
    REAL* vertices; // Polygon must be defined in an oriented fashion. See https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
    REAL maxRes;
    REAL minRes;
    REAL decayRes;
    bool indicator;

    PolygonTest(int nVertices, REAL* vertices, REAL res) :
    nVertices(nVertices),
    vertices(vertices),
    maxRes(1.0), // unused
    minRes(res),
    decayRes(1.0), // unused
    indicator(true)
    {
        assert(nVertices >= 3);
    }

    PolygonTest(int nVertices, REAL* vertices, REAL maxRes, REAL minRes, REAL decayRes) :
    nVertices(nVertices),
    vertices(vertices),
    maxRes(maxRes),
    minRes(minRes),
    decayRes(decayRes),
    indicator(false)
    {
        assert(nVertices >= 3);
    }

    REAL operator()(REAL x, REAL y) {
        // Check if point is inside polygon
        int j = nVertices - 1;
        bool inside = false;
        bool flippedYet = false;

        for (int i = 0; i < nVertices; i++) {
            if ((vertices[2 * i + 1] > y) != (vertices[2 * j + 1] > y)) {
                REAL xx = (vertices[2 * j] - vertices[2 * i]) * (y - vertices[2 * i + 1]) / (vertices[2 * j + 1] - vertices[2 * i + 1]) + vertices[2 * i];
                if (x < xx) {
                    inside = !inside;
                }
            }
            j = i;
        }

        if (inside) {
            return minRes;
        }

        // Otherwise, return outside resolution
        if (indicator) {
            return INFINITY;
        }
        else {
            REAL res = maxRes - (maxRes - minRes) * exp(-distanceToSegment(x, y, vertices[2 * (nVertices - 1)], vertices[2 * (nVertices - 1) + 1], vertices[0], vertices[1]) / decayRes);
            for (int i = 1; i < nVertices; i++) {
                REAL testres = maxRes - (maxRes - minRes) * exp(-distanceToSegment(x, y, vertices[2 * (i - 1)], vertices[2 * (i - 1) + 1], vertices[2 * i], vertices[2 * i + 1]) / decayRes);
                if (testres < res) { res = testres; }
            }

            return res;
        }
    }
};


// By default, return a maximum "resolution" that will always be large.
REAL defaultFunction(REAL x, REAL y, REAL* args) {
    return INFINITY;
}


struct Callable {
    std::vector<REAL*> arg_list;
    std::vector<REAL (*)(REAL, REAL, REAL*)> function_list;
    std::vector<PointTest> point_list;
    std::vector<LineTest> line_list;
    std::vector<PolygonTest> polygon_list;

    REAL* args; // default arguments, never used

    REAL operator()(REAL x, REAL y) {
        // Return the default function if there are no constraints defined
        if (function_list.size() == 0 && point_list.size() == 0 && line_list.size() == 0 && polygon_list.size() == 0) {
            return defaultFunction(x, y, args);
        }

        REAL maxres = INFINITY;

        // Loop through the defined functions and take their minimum envelope
        for (int i = 0; i < function_list.size(); i++) {
            REAL testres = function_list[i](x, y, arg_list[i]);
            if (testres < maxres) { maxres = testres; }
        }

        // Loop through the points
        for (int i = 0; i < point_list.size(); i++) {
            REAL testres = point_list[i](x, y);
            if (testres < maxres) { maxres = testres; }
        }

        // Loop through the lines
        for (int i = 0; i < line_list.size(); i++) {
            REAL testres = line_list[i](x, y);
            if (testres < maxres) { maxres = testres; }
        }

        // Loop through the polygons
        for (int i = 0; i < polygon_list.size(); i++) {
            REAL testres = polygon_list[i](x, y);
            if (testres < maxres) { maxres = testres; }
        }

        return maxres;
    }
} callable;


extern "C" int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area)
{
    REAL dxoa, dxda, dxod;
    REAL dyoa, dyda, dyod;
    REAL oalen, dalen, odlen;
    REAL maxlen;

    dxoa = triorg[0] - triapex[0];
    dyoa = triorg[1] - triapex[1];
    dxda = tridest[0] - triapex[0];
    dyda = tridest[1] - triapex[1];
    dxod = triorg[0] - tridest[0];
    dyod = triorg[1] - tridest[1];
    /* Find the squares of the lengths of the triangle's three edges. */
    oalen = dxoa * dxoa + dyoa * dyoa;
    dalen = dxda * dxda + dyda * dyda;
    odlen = dxod * dxod + dyod * dyod;
    /* Find the square of the length of the longest edge. */
    maxlen = (dalen > oalen) ? dalen : oalen;
    maxlen = (odlen > maxlen) ? odlen : maxlen;

    REAL trix = (triorg[0] + tridest[0] + triapex[0]) / 3.0;
    REAL triy = (triorg[1] + tridest[1] + triapex[1]) / 3.0;
    REAL maxres = callable(trix, triy);

    if (maxlen > maxres * maxres) { return 1; }
    else { return 0; }
}


extern "C" void setTriangleResolutionFunction(REAL (*f)(REAL, REAL, REAL*), REAL* args) {
    callable.function_list.push_back(f);
    callable.arg_list.push_back(args);
    return;
}


extern "C" void setTriangleResolutionPointOffsetFunction(REAL maxRes, REAL minRes, REAL decayRes, REAL x0, REAL y0, REAL offset) {
    PointTest test(maxRes, minRes, decayRes, x0, y0, offset);
    callable.point_list.push_back(test);
    return;
}


extern "C" void setTriangleResolutionPointFunction(REAL maxRes, REAL minRes, REAL decayRes, REAL x0, REAL y0) {
    PointTest test(maxRes, minRes, decayRes, x0, y0);
    callable.point_list.push_back(test);
    return;
}


extern "C" void setTriangleResolutionLineOffsetFunction(REAL maxRes, REAL minRes, REAL decayRes, REAL xi, REAL yi, REAL xf, REAL yf, REAL offset) {
    LineTest test(maxRes, minRes, decayRes, xi, yi, xf, yf, offset);
    callable.line_list.push_back(test);
    return;
}


extern "C" void setTriangleResolutionLineFunction(REAL maxRes, REAL minRes, REAL decayRes, REAL xi, REAL yi, REAL xf, REAL yf) {
    LineTest test(maxRes, minRes, decayRes, xi, yi, xf, yf);
    callable.line_list.push_back(test);
    return;
}


extern "C" void setTriangleResolutionPolygonFunction(REAL maxRes, REAL minRes, REAL decayRes, int nVertices, REAL* vertices) {
    PolygonTest test(nVertices, vertices, maxRes, minRes, decayRes);
    callable.polygon_list.push_back(test);
    return;
}


extern "C" void setTriangleResolutionPolygonIndicatorFunction(int nVertices, REAL* vertices, REAL res) {
    PolygonTest test(nVertices, vertices, res);
    callable.polygon_list.push_back(test);
    return;
}


extern "C" void clearTriangleResolutionFunction() {
    callable.function_list.clear();
    callable.arg_list.clear();
    callable.polygon_list.clear();
    callable.line_list.clear();
    callable.point_list.clear();
    return;
}

template <typename T>
class nacl{

public:
    const T boxlength;
    const int n;
    const T x[8];
    const T y[8];
    const T z[8];
    const T q[8];

    nacl() : boxlength(2.0), n(8),
        x
        {
        0.5,
        1.5,
        0.5,
        1.5,
        0.5,
        1.5,
        0.5,
        1.5
        },
        y
        {
        0.5,
        0.5,
        1.5,
        1.5,
        0.5,
        0.5,
        1.5,
        1.5
        },
        z
        {
        0.5,
        0.5,
        0.5,
        0.5,
        1.5,
        1.5,
        1.5,
        1.5
        },
        q
        {
        1.0,
        -1.0,
        -1.0,
        1.0,
        -1.0,
        1.0,
        1.0,
        -1.0
        }
    {}
};


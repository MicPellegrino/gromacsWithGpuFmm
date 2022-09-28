template <typename T>
class random_parts_2{

public:
    const T boxlength;
    const int n;
    const T x[2];
    const T y[2];
    const T z[2];
    const T q[2];

    random_parts_2() : boxlength(1.0),n(2),
    x
    {
    0.09375,
    0.90625
    },
    y
    {
    0.5,
    0.5
    },
    z
    {
    0.5,
    0.5
    },
    q
    {
    -1.0,
    1.0
    }
    {}
};

template <typename T>
class parts_2{

public:
    const T boxlength;
    const int n;
     T x[2];
     T y[2];
     T z[2];
     T q[2];

    parts_2() : boxlength(2.0),n(2),
    x
    {
    0.5,
    1.5
    },
    y
    {
    0.5,
    0.5
    },
    z
    {
    0.5,
    0.5
    },
    q
    {
    -1.0,
    1.0
    }
    {}
};




#pragma once

union v2_x4
{
    struct
    {
        v1_x4 x, y;
    };

    v1_x4 e[2];
};


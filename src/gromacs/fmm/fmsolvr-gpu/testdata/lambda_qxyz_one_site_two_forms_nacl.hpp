const REAL boxlength=2;
const int num_of_sites=1;
const int site_ids[2]=
{
1,
2
};
const int num_of_forms_per_site[2]=
{
1,
2
};
const int n_0_0 = 0;
const REAL l_0_0 = 1.0;
const REAL x_0_0[n_0_0] = 
{
};
const REAL y_0_0[n_0_0] = 
{
};
const REAL z_0_0[n_0_0] = 
{
};
const REAL q_0_0[n_0_0] = 
{
};
const REAL **x_site0 = new const REAL *[1];
const REAL **y_site0 = new const REAL *[1];
const REAL **z_site0 = new const REAL *[1];
const REAL **q_site0 = new const REAL *[1];
const int n_1_0 = 8;
const REAL l_1_0 = 0.1;
const REAL x_1_0[n_1_0] = 
{
0.5,
1.5,
0.5,
1.5,
0.5,
1.5,
0.5,
1.5
};
const REAL y_1_0[n_1_0] = 
{
0.5,
0.5,
1.5,
1.5,
0.5,
0.5,
1.5,
1.5
};
const REAL z_1_0[n_1_0] = 
{
0.5,
0.5,
0.5,
0.5,
1.5,
1.5,
1.5,
1.5
};
const REAL q_1_0[n_1_0] = 
{
1,
-1,
-1,
1,
-1,
1,
1,
-1
};
const int n_1_1 = 8;
const REAL l_1_1 = 0.9;
const REAL x_1_1[n_1_1] = 
{
0.5,
1.5,
0.5,
1.5,
0.5,
1.5,
0.5,
1.5
};
const REAL y_1_1[n_1_1] = 
{
0.5,
0.5,
1.5,
1.5,
0.5,
0.5,
1.5,
1.5
};
const REAL z_1_1[n_1_1] = 
{
0.5,
0.5,
0.5,
0.5,
1.5,
1.5,
1.5,
1.5
};
const REAL q_1_1[n_1_1] = 
{
1,
-1,
-1,
1,
-1,
1,
1,
-1
};
const REAL **x_site1 = new const REAL *[2];
const REAL **y_site1 = new const REAL *[2];
const REAL **z_site1 = new const REAL *[2];
const REAL **q_site1 = new const REAL *[2];

const REAL ***x_data = new const REAL **[1];
const REAL ***y_data = new const REAL **[1];
const REAL ***z_data = new const REAL **[1];
const REAL ***q_data = new const REAL **[1];
int **n_data = new int *[1];
int *n_site0 = new int[1];
int *n_site1 = new int[2];

REAL **l_data = new REAL *[1];
REAL *l_site0 = new REAL[1];
REAL *l_site1 = new REAL[2];


int acu_data(){

x_site0[0] = x_0_0;
y_site0[0] = y_0_0;
z_site0[0] = z_0_0;
q_site0[0] = q_0_0;
n_site0[0] = n_0_0;
l_site0[0] = l_0_0;
x_site1[0] = x_1_0;
y_site1[0] = y_1_0;
z_site1[0] = z_1_0;
q_site1[0] = q_1_0;
n_site1[0] = n_1_0;
l_site1[0] = l_1_0;
x_site1[1] = x_1_1;
y_site1[1] = y_1_1;
z_site1[1] = z_1_1;
q_site1[1] = q_1_1;
n_site1[1] = n_1_1;
l_site1[1] = l_1_1;
x_data[0] = x_site0;
y_data[0] = y_site0;
z_data[0] = z_site0;
q_data[0] = q_site0;
n_data[0] = n_site0;
l_data[0] = l_site0;
x_data[1] = x_site1;
y_data[1] = y_site1;
z_data[1] = z_site1;
q_data[1] = q_site1;
n_data[1] = n_site1;
l_data[1] = l_site1;

return 1;
}

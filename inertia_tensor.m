function I = inertia_tensor(varargin)
% This function is an implementation of the "List of 3D inertia tensors" in
% wikipedia page https://en.wikipedia.org/wiki/List_of_moments_of_inertia.
%
% Usage:
%   I = inertia_tensor(input_arguments)
%
%   Input arguments format:
%       [name_of_the_shape] + [parameters] + <optional parameters>
%
% Example:
%   I = inertia_tensor('solid_sphere', 5)
%   I = inertia_tensor('cub', 3, 4, 5)
%   I = inertia_tensor('cyl', 1, 5, 'mass', 10)
%
%   -------------------------------------------------------------------
%          Name of the shape      |              parameters
%   -------------------------------------------------------------------
%             Solid_Sphere        |              radius (r)
%            Hollow_Sphere        |              radius (r)
%            Solid_Ellipsoid      |          semi-axes (a,b,c)
%          Right_Circular_Cone    |        radius (r), height (h)
%             solid_CUBoid        |  width (w), height (h), depth (d)
%           Slender_Rod_End       |              length (l)
%           Slender_Rod_Cen       |              length (l)
%            solid_CYLinder       |       radius (r), height (h)
%           Thick_Walled_Tube     | inner/outer radius (r1, r2), height (h)

if (nargin <= 1)
    error("Too few input arguments. Please refer to the help of the function.");
end

shape_type = varargin{1};
if ~ischar(shape_type)
    error("The first argument does not specify the name of the shape.");
end

% Define the mass of the shape
mass = 1.0;
for i=2:length(varargin)
    if (strcmpi(varargin{i}, 'mass')) && isnumeric(varargin{i+1})
        mass = varargin{i+1};
        break
    end    
end

switch lower(shape_type)
    case {'solid_sphere', 'ss'}
        r = varargin{2};
        I = inertia_ss(r, mass);
    case {'hollow_sphere', 'hs'}
        r = varargin{2};
        I = inertia_hs(r, mass);
    case {'solid_ellipsoid', 'se'}
        a = varargin{2};
        b = varargin{3};
        c = varargin{4};
        I = inertia_se(a, b, c, mass);
    case {'right_circular_cone', 'rcc'}
        r = varargin{2};
        h = varargin{3};
        I = inertia_rcc(r, h, mass);
    case {'solid_cuboid', 'cub'}
        w = varargin{2};
        h = varargin{3};
        d = varargin{4};
        I = inertia_cub(w, h, d, mass);
    case {'slender_rod_end', 'sre'}
        l = varargin{2};
        I = inertia_sre(l, mass);
    case {'slender_rod_cen', 'src'}
        l = varargin{2};
        I = inertia_src(l, mass);
    case {'solid_cylinder', 'cyl'}
        r = varargin{2};
        h = varargin{3};
        I = inertia_cyl(r, h, mass);
    case {'thick_walled_tube', 'twt'}
        r1 = varargin{2};
        r2 = varargin{3};
        h = varargin{4};
        I = inertia_twt(r1, r2, h, mass);
end


function I = inertia_ss(r, m)
I = [2/5*m*r*r, 0, 0; 0, 2/5*m*r*r, 0; 0, 0, 2/5*m*r*r];

function I = inertia_hs(r, m)
I = [2/3*m*r*r, 0, 0; 0, 2/3*m*r*r, 0; 0, 0, 2/3*m*r*r];

function I = inertia_se(a, b, c, m)
I = [1/5*m*(b*b+c*c), 0, 0; 0, 1/5*m*(a*a+c*c), 0; 0, 0, 1/5*m*(a*a+b*b)];

function I = inertia_rcc(r, h, m)
I = [3/5*m*h*h+3/20*m*r*r, 0, 0; 0, 3/5*m*h*h+3/20*m*r*r, 0; 0, 0, 3/10*m*r*r];

function I = inertia_cub(w, h, d, m)
I = [1/12*m*(h*h+d*d), 0, 0; 0, 1/12*m*(w*w+d*d), 0; 0, 0, 1/12*m*(w*w+h*h)];

function I = inertia_sre(l, m)
I = [1/3*m*l*l, 0, 0; 0, 0, 0; 0, 0, 1/3*m*l*l];

function I = inertia_src(l, m)
I = [1/12*m*l*l, 0, 0; 0, 0, 0; 0, 0, 1/12*m*l*l];

function I = inertia_cyl(r, h, m)
I = [1/12*m*(3*r*r+h*h), 0, 0; 0, 1/12*m*(3*r*r+h*h), 0; 0, 0, 1/2*m*r*r];

function I = inertia_twt(r1, r2, h, m)
I = [1/12*m*(3*(r1*r1+r2*r2)+h*h), 0, 0; 0, 1/12*m*(3*(r1*r1+r2*r2)+h*h), 0; 0, 0, 1/2*m*(r1*r1+r2*r2)];
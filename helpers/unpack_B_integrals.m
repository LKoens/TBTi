function out = unpack_B_integrals(res)
%% UNPACK_B_INTEGRALS will take the vector output used in the construction of
% B and assign its components to a struct with correct field names.
out = struct();

% Stokeslet.
out.Bi  = res(1);
out.B11 = res(2);
out.B12 = res(3);
out.B13 = res(4);
out.B22 = res(5);
out.B23 = res(6);
out.B33 = res(7);

% Stokeslet boundary correction.
out.Biw  = res(8);
out.B11w = res(9);
out.B12w = res(10);
out.B13w = res(11);
out.B22w = res(12);
out.B23w = res(13);
out.B33w = res(14);

% Source dipole.
out.Biwsd  = res(15);
out.B11wsd = res(16);
out.B12wsd = res(17);
out.B13wsd = res(18);
out.B22wsd = res(19);
out.B23wsd = res(20);
out.B33wsd = res(21);

% Stresslet.
out.B13wss = res(22);
out.B23wss = res(23);

end
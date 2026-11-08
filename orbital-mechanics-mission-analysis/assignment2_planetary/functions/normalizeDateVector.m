function dateN = normalizeDateVector(dateIn)
%NORMALIZEDATEVECTOR Normalizes a date vector in the format [Y M D h m s].
%
% PROTOTYPE:
%   dateN = normalizeDateVector(dateIn)
%
% PURPOSE:
%   Brings seconds, minutes, and hours into their standard ranges:
%     - s in [0, 59]
%     - m in [0, 59]
%     - h in [0, 23]
%   by carrying the excess to the next higher unit.
%   This is useful, for example, to handle inputs like "... 05:24:60".
%
% INPUT:
%   dateIn [1x6]  Date vector: [Y M D h m s]
%
% OUTPUT:
%   dateN  [1x6]  Normalized date vector: [Y M D h m s]
%
% NOTE:
%   The normalization handles ONLY carries s->m, m->h, h->D.
%   It does not handle month/year changes (D beyond end of month), because in the project
%   it is not needed for the cases considered.
%
% AUTHORS:
%   Simone Frescura
%   Nicola Pifano
%   Luca Puzzi
%   Lorenzo Veronese
%==========================================================================

Y = dateIn(1);  M = dateIn(2);  D = dateIn(3);
h = dateIn(4);  m = dateIn(5);  s = dateIn(6);

% Seconds -> minutes carry
carry_m = floor(s/60);
s = s - 60*carry_m;
m = m + carry_m;

% Minutes -> hours carry
carry_h = floor(m/60);
m = m - 60*carry_h;
h = h + carry_h;

% Hours -> days carry
carry_d = floor(h/24);
h = h - 24*carry_d;
D = D + carry_d;

dateN = [Y M D h m s];

end


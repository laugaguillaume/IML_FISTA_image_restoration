function make(flag, outdir)
%
% Usage:
%    make         ->  compile all
%    make openmp  ->  compile all with OpenMP (for compilers with OpenMP support)
%
% HINT: 
%  remember to add the folder "build\cpu" or "build\openmp" in MATLAB path


% basic
include = '"-Iinclude/cpp -Iinclude/Eigen"';

% default
if nargin < 1
    flag = '';
end
if nargin < 2
    outdir = 'build';
end

% read the inputs
if strcmpi(flag, 'openmp')
    outdir = [outdir '/openmp'];
    opts = ['-outdir ' outdir];
    opts = [opts ' COMPFLAGS="$COMPFLAGS /openmp"'];
else
    outdir = [outdir '/cpu'];
    opts = ['-outdir ' outdir];
end
opts = ['"' opts '"'];

% make the folder
mkdir(outdir);



% find the the directories
paths = genpath('NLST_algo/conv_opt');
dirs = regexp(paths, pathsep, 'split');

% look for the files to be compiled
files = discover_files(dirs,'cpp');

% compile
for n=1:length(files)
    
    % extract the name
    str = regexp(files{n},['([^\' filesep ']+$)'], 'match');   % {<name>.cpp}
    str = regexp(str{1},'\w+','match');                        % {<name> cpp}
    name = str{1};

    fprintf('Compiling %s... ', [name '.' str{2}]);
    
    % compile
    n
    pwd
    files{n}
    mex(files{n}, ['"-output ' name '"'], opts, include);
    
    fprintf('done.\n');
end





 function images = discover_files( dir_path, suffix, N_max )
%function images = discover_files( dir_path, suffix, N_max )
%
%  Created on: 09/06/10 - Giovanni Chierchia
% Revision #1: 05/09/10 - Giovanni Chierchia (BUGFIX: limite max dei file)
%
%
% La funzione restituisce tutti file che terminano con la stringa "suffix" 
% presenti nelle cartelle "dir_path{1:end}", . Il parametro opzionale 
% "N_max" impone il limite massimo sul numero totale di file selezionati.



% controllo dei parametri
if ~ischar( suffix ) || length(suffix) < 3
    error('il parametro "suffix" deve essere una stringa di almeno 3 caratteri');
end
if ischar( dir_path ) 
    path = { dir_path };
elseif iscell( dir_path )
    path = dir_path;
else
    error('il parametro "dir_path" deve essere una stringa o una cella di stringhe');
end

% parametri di default
if nargin < 3
    N_max = 0;
end

% immagini per singola directory
N_path = length(path);
if N_max ~= 0
    N_limit = ceil(N_max / N_path);
end

% scorri le directory
idx = 1;
for j=1:N_path

    if ~ischar( path{j} )
        error('path non corretto');
    end
    
    % elenca tutti i file nella directory
    files = dir( path{j} );

    % scorri i file
    count = 0;
    for i=1:length(files)                          
    
        % nome del file
        filename = files(i).name;
    
        % prendi il file se ha l'estensione giusta
        sfx = length(suffix);
        if length(filename) >= sfx && strcmpi( filename(end-sfx+1:end), suffix )
            images{idx} = fullfile( path{j}, filename );
            idx = idx+1;
            count = count + 1;
        end
    
        % non prendere più di 'N_limit' immagini per ogni directory
        if N_max ~= 0 && (count == N_limit || idx-1 == N_max)
            break;
        end
    
    end
    
    % non prendere più di 'N_max' immagini
    if N_max ~= 0 && idx-1 == N_max
        break;
    end
    
end

% assegna l'uscita se non sono state trovate immagini
if idx == 1
    images = {};
end
classdef SPASM < handle
% SPASM Builds an app with the ability to detect binding events from noisy single
% molecule data and produce relevant figures including ensemble averages.
%   SPASM creates a new app instance and displays a GUI for interacting with the
%   app.
%
%   SPASM(txtInputFcn) accepts either a function handle or a character vector
%   containing the filepath of a MATLAB function. This function will be used to read .txt
%   files. Refer to the user guide for more information.
%
%   SPASM(txtInputFcn, txtOutputFcn) also accepts either a function handle or a
%   character vector containing the filepath of a MATLAB function which will be used to
%   print data to .txt files. Refer to the user guide for more information.
%
%   myApp = SPASM( ___ ) returns the app instance in myApp. Use dot indexing to
%   access properties and methods of myApp.
% 
%   See also properties, methods.
%
% A tutorial can be found in the user guide, which is located at:
% https://github.com/GreenbergLab/SPASM
%
% This program uses the Signal Processing Toolbox to smooth the covariance with a
% Savitzky-Golay filter. The Optimization Toolbox is used to fit exponential curves to the
% ensemble averages. If either toolbox is not installed, the program will run without the
% corresponding functionality.
%
% This program should be compatible with MATLAB releases R2017b through at least R2020a.

    % These are observable - a function is called when their values change.
    properties (SetObservable, AbortSet)
        allEvents            % A double array containing the start (last pt before event) and end (last pt of event) times for all events (both active and removed).
        activeEvents         % A logical array that indicates which events are active (i.e. not removed automatically by program or manually by user).
        minSep = 173         % Events separated by less than minSep number of points are removed automatically.
        minDur = 150         % Events shorter than minDur number of points are removed automatically.
        autop1               % The lower peak of the covariance histogram, which is used in the peak-to-peak method.
        autop2               % The upper peak of the covariance histogram, which is used in the peak-to-peak method.
        automin              % The minimum value between the two peaks of the covariance histogram, which is used in the threshold method.
        manp1                % The user has the option to specify peak 1 manually, in case the program does not choose an appropriate value.
        manp2                % The user has the option to specify peak 2 manually, in case the program does not choose an appropriate value.
        manmin               % The user has the option to specify the minimum manually, in case the program does not choose an appropriate value.
        useman = false       % A boolean determining if the program will use the manual peaks and minimum (as opposed to the automatic values).
        usemin = false       % A boolean determining if the program will use the threshold method (as opposed to the peak-to-peak method).
        flipped = false      % A boolean determining if the raw data is upside down and needs to be flipped.
        showEvents = true    % A boolean determining if red patches on the main axes which indicate events should be shown.
        showRemoved = false  % A boolean determining if the user should be able to view removed events on the individual event/likelihood axes, as well as the patches for removed events on the main axes when the user mouses over them.
    end
    
    % Remaining properties are contained in these structs.
    properties
        load      % Contains properties related to the load panel in tab 1.
        analyze   % Contains properties related to the analyze panel in tab 1.
        controls  % Contains properties related to the controls panel in tab 1 (Find Events button, e.g.).
        excel     % Contains properties related to the save to excel panel in tab 1.
        misc      % Contains miscellaneous properties (the main figure, e.g.).
        data      % Contains relevant data for tab 1 (bead positions, sampling frequency, ensemble averages, e.g.).
        load2     % Contains properties related to the load panel in tab 2.
        analyze2  % Contains properties related to the analyze panel in tab 2.
        data2     % Contains relevant data for tab 2 (event durations, number of events, e.g.).
    end

    % Public non-static methods
    methods
        
        function app = SPASM(txtInputFcn, txtOutputFcn)
            % Constructor.
            
            try
                app.initialize();
                app.misc.initBackupTab1 = app.storeBackupTab1();
                app.misc.initBackupTab2 = app.storeBackupTab2();
                
                if nargin >= 1 && ~isdeployed
                    % If user supplied a custom input fcn:
                    if isa(txtInputFcn, 'function_handle')
                        % If they supplied a function handle, store the handle.
                        app.load.txtInputFcn = txtInputFcn;
                    else
                        % If they did not supply a function handle, assume it is a character
                        % vector containing the name of a .m file:
                        try
                            [path, name] = fileparts(txtInputFcn);
                            oldDir = pwd; % Store present directory.
                            if ~isempty(path)
                                cd(path); % Change directory to wherever the .m file is located.
                            end
                            app.load.txtInputFcn = str2func(name); % Obtain a function handle for the supplied .m file.
                            cd(oldDir); % Return to original directory.
                        catch
                            % An error could arise for any number of reasons, including if the
                            % file is not found. Use the default function and warn the user.
                            warndlg(sprintf('Unable to convert first input argument to a function. Using the default data input function to read .txt files.'))
                            app.load.txtInputFcn = @app.defaultTxtInput;
                        end
                    end
                    try
                        % Try to run the function...
                        app.load.txtInputFcn();
                    catch ME
                        % ...catch the inevitable error...
                        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                            % ...and, if the error indicates that the function is not defined,
                            % just use the default function and warn the user.
                            warndlg(sprintf('Unable to find %s(). Using the default data input function to read .txt files.', char(txtInputFcn)))
                            app.load.txtInputFcn = @app.defaultTxtInput;
                        end
                    end
                end

                if nargin >= 2 && ~isdeployed
                    % If user supplied a custom output fcn:
                    if isa(txtOutputFcn, 'function_handle')
                        % If they supplied a function handle, store the handle.
                        app.load.txtOutputFcn = txtOutputFcn;
                    else
                        % If they did not supply a function handle, assume it is a character
                        % vector containing the name of a .m file:
                        try
                            [path, name] = fileparts(txtOutputFcn);
                            oldDir = pwd; % Store present directory.
                            if ~isempty(path)
                                cd(path); % Change directory to wherever the .m file is located.
                            end
                            app.load.txtOutputFcn = str2func(name); % Obtain a function handle for the supplied .m file.
                            cd(oldDir); % Return to original directory.
                        catch
                            % An error could arise for any number of reasons, including if the
                            % file is not found. Use the default function and warn the user.
                            warndlg(sprintf('Unable to convert second input argument to a function. Using the default txt output function to write .txt files.'))
                            app.load.txtOutputFcn = @app.defaultTxtOutput;
                        end
                    end
                    try
                        % Try to run the function...
                        app.load.txtOutputFcn();
                    catch ME
                        % ...catch the inevitable error...
                        if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                            % ...and, if the error indicates that the function is not defined,
                            % just use the default function and warn the user.
                            warndlg(sprintf('Unable to find %s(). Using the default txt output function to write .txt files.', char(txtOutputFcn)))
                            app.load.txtOutputFcn = @app.defaultTxtOutput;
                        end
                    end
                end

                % If user does not indicate that they want to save the app, clear it to avoid
                % any output.
                if nargout == 0
                    clear app
                end
            catch ME
                app.unexpectedError(ME, sprintf('An unexpected error occurred while trying to start the app:\n\n%s', ME.message));
                try
                    close(app.misc.fig);
                catch
                end
                delete(app);
                
                % If user does not indicate that they want to save the app, clear it to
                % avoid any output.
                if nargout == 0
                    clear app
                end
            end
        end
        
        function delete(app)
            % Destructor.
            
            % The delete method of the handle superclass will delete the object, but we
            % want to make sure no figure is left behind once the object is deleted. If a
            % figure were left behind, buttons e.g. would be clickable but their callbacks
            % would throw errors. Deleting the figure will close it and delete all of its
            % children (buttons, axes, etc).
            try
                delete(app.misc.fig);
            catch
            end
        end
        
        function trimBeads = getTrimBeads(app)
            % GETTRIMBEADS Returns the current beads.
            %   trimBeads = GETTRIMBEADS(app) returns app's current beads, excluding data
            %   which was removed by the user.
            
            trimBeads = app.data.origBeads(app.controls.deselect.indices, :);
        end
        
        function trimData = getTrimData(app)
            % GETTRIMDATA Returns the current data table.
            %   trimData = GETTRIMDATA(app) returns app's current data table, excluding
            %   rows that correspond to data which was removed by the user. This table is
            %   the same table returned by app's .txt input function.
            %
            % See also defaulttxtinput.
            
            trimData = app.data.origData(app.controls.deselect.indices, :);
        end
        
        function [A, B] = getActiveBeads(app)
            % GETACTIVEBEADS Returns the position of each bead separately.
            %   [A, B] = GETACTIVEBEADS(app) returns the positions of app's current beads
            %   separately, excluding data which was removed by the user.
            
            trimBeads = app.getTrimBeads;
            if ~isempty(trimBeads)
                A = trimBeads(:,1);
                B = trimBeads(:,2);
            else
                A = [];
                B = [];
            end
        end
        
        function avg = getBeadsToAnalyze(app)
            % GETBEADSTOANALYZE Returns the average of the two beads' positions.
            %   avg = GETBEADSTOANALYZE(app) returns the averaged position of app's two
            %   beads, excluding data which was removed by the user. If the user has
            %   specified that only one bead should be analyzed, that bead's position will
            %   be returned instead.
            
            [A, B] = app.getActiveBeads;
            if app.controls.beadABtn.Value
                % Just bead A.
                avg = A;
            elseif app.controls.beadBBtn.Value
                % Just bead B.
                avg = B;
            else
                % Both.
                avg = (A+B)/2;
            end
        end
        
        function events = getActiveEvents(app)
            % GETACTIVEEVENTS Returns the currently active events.
            %   events = GETACTIVEEVENTS(app) returns app's currently active events in a
            %   2-column array: the first column gives the event starts (the last index
            %   after each event within the app's beads), the second column gives the
            %   event ends (the last index of each event), and the rows are ordered
            %   chronologically.
            
            if ~isempty(app.allEvents) && ~isempty(app.activeEvents)
                events = app.allEvents(app.activeEvents,:);
            else
                events = [];
            end
        end
        
        function dur = getActiveDur(app)
            % GETACTIVEDUR Returns the durations of the currently active events.
            %   dur = GETACTIVEDUR(app) returns the durations, in number of points, of
            %   app's currently active events. dur(i) gives the duration of event i, with
            %   the events ordered chronologically.
            %
            % You can access the sampling frequency with app.data.fs.
            
            events = app.getActiveEvents;
            if ~isempty(events)
                dur = events(:,2) - events(:,1);
            else
                dur = [];
            end
        end
        
        function sep = getActiveSep(app)
            % GETACTIVESEP Returns the number of points between consecutive events for the
            % currently active events.
            %   sep = GETACTIVESEP(app) returns the number of points separating app's
            %   currently active consecutive events. sep(i,1) gives the number of points
            %   before event i, and sep(i,2) gives the number of points after event i,
            %   with the events ordered chronologically.
            %
            % You can access the sampling frequency with app.data.fs.
            
            events = app.getActiveEvents;
            if ~isempty(events)
                sep = [events(:,1); length(app.getTrimBeads)] - [1; events(:,2)];
            else
                sep = [];
            end
        end
        
        function [pos, steps] = getPosNStepSizes(app, whichBead)
            % GETPOSNSTEPSIZES Returns the average position of a bead just before and just
            % after binding and unbinding, as well as the estimated step sizes.
            %   [pos, steps] = GETPOSNSTEPSIZES(app, whichBead) returns the average
            %   position, in nm, of the requested bead just before and just after each
            %   currently active event's onset and conclusion, as well as the estimated
            %   step sizes based on these positions. If whichBead equals 'B', this
            %   function will analyze bead B. Otherwise, this function will analyze bead
            %   A.
            %
            %   pos is a 4-column array, where row i corresponds to event i, with the
            %   events ordered chronologically:
            %     The first column contains the position of the bead averaged over 10 ms
            %       just prior to each event's onset.
            %     The second column contains the position of the bead averaged over 10 ms
            %       just after each event's onset.
            %     The third and fourth columns contain similar values for each event's
            %       conclusion.
            %
            %   steps is a 3-column array, where row i corresponds to event i, with the
            %   events ordered chronologically:
            %     The first column contains estimates for the size of the requested bead's
            %       first substep during each event.
            %     The second column contains estimates for the bead's second substep
            %       during each event.
            %     The third column contains estimates for the total step (substep 1 plus
            %       substep 2).
            
            events = app.getActiveEvents;
            [A, B] = app.getActiveBeads;
            if app.flipped
                A = -A;
                B = -B;
            end
            if strcmp(whichBead, 'B')
                bead = B;
            else
                bead = A;
            end
            N = size(events, 1);
            M = numel(bead);
            posBefB = zeros(N, 1);
            posAftB = zeros(N, 1);
            posBefU = zeros(N, 1);
            posAftU = zeros(N, 1);
            events = [0 0; events; M, M];
            for i = 2:N+1
                posBefB(i-1) = mean(bead(max([events(i-1,2)+1, events(i,1)-app.data.fs/100]):events(i,1)));
                posAftB(i-1) = mean(bead(events(i,1)+1:min([events(i,1)+1+app.data.fs/100, events(i,2)])));
                posBefU(i-1) = mean(bead(max([events(i,2)-app.data.fs/100, events(i,1)+1]):events(i,2)));
                posAftU(i-1) = mean(bead(events(i,2)+1:min([events(i,2)+1+app.data.fs/100, events(i+1,1)])));
            end
            step1 = posAftB - posBefB;
            totalstep = posBefU - posAftU;
            step2 = totalstep - step1;
            pos = [posBefB...
                posAftB...
                posBefU...
                posAftU];
            steps = [step1...
                step2...
                totalstep];
        end
        
        function force = getForce(app, whichBead)
            % GETFORCE Returns the force on the requested bead before, during, and after
            % each event.
            %   force = GETFORCE(app, whichBead) returns the force, in pN, of the
            %   requested bead before, during, and after each currently active event. If
            %   whichBead equals 'B', this function will analyze bead B. Otherwise, this
            %   function will analyze bead A.
            %
            %   force is a 5-column array, where row i corresponds to event i, with the
            %   events ordered chronologically:
            %     The first column contains the average force on the bead over the entire
            %       duration separating a given event from the preceding event (or, in the
            %       case of the first event, from the beginning of the file).
            %     The second column contains the average force on the bead during the
            %       event.
            %     The third column contains the average force on the bead over the entire
            %       duration separating a given event from the following event (or, in the
            %       case of the final event, from the end of the file).
            %     The fourth and fifth columns contain the force on the bead averaged over
            %       10 ms just before and just after each event's conclusion.
            
            events = app.getActiveEvents;
            sep = app.getActiveSep;
            [A, B] = app.getActiveBeads;
            if app.flipped
                A = -A;
                B = -B;
            end
            if strcmp(whichBead, 'B')
                bead = B;
                K = app.load.KB;
            else
                bead = A;
                K = app.load.KA;
            end
            N = size(events, 1);
            M = numel(bead);
            FbefE = zeros(N, 1);
            FdurE = zeros(N, 1);
            FaftE = zeros(N, 1);
            FbefD = zeros(N, 1);
            FaftD = zeros(N, 1);
            events = [events; M, M];
            for i = 1:N
                FbefE(i) = mean(bead(events(i,1)+1-sep(i):events(i,1))*K);
                FdurE(i) = mean(bead(events(i,1)+1:events(i,2))*K);
                FaftE(i) = mean(bead(events(i,2)+1:events(i,2)+sep(i+1))*K);
                FbefD(i) = mean(bead(max([1, events(i,2)-app.data.fs/100, events(i,1)+1]):events(i,2))*K);
                FaftD(i) = mean(bead(events(i,2)+1:min([events(i,2)+1+app.data.fs/100, events(i+1,1)]))*K);
            end
            force = [FbefE...
                FdurE...
                FbefD...
                FaftD...
                FaftE];
        end
    end
    
    % Public static methods (these do not require an instance of this class)
    methods (Static)
        
        function [fs, KA, KB, CALA, CALB, tData, header] = defaultTxtInput(file)
            % DEFAULTTXTINPUT Reads .txt file and returns the contents.
            %   [fs, KA, KB, CALA, CALB, tData, header] = DEFAULTTXTINPUT(file) assumes
            %   file is a text file with the following format:
            %
            %       key1        value1
            %       ...         ...
            %       keyN        valueN
            %       colName1    colName2    ...     colNameM
            %       data(1,1)   data(1,2)   ...     data(1,M)
            %       ...         ...         ...     ...
            %       data(L,1)   data(L,2)   ...     data(L,M)
            %
            %   where columns of text are tab-delimited. The following must be true:
            %
            %   (1) All data must be numeric text (e.g. '1.0'), and the first row of data
            %   must be the first row in the file which contains only numeric text.
            %
            %   (2) There must be keys labeled 'K1', 'K3', 'CAL1', 'CAL3', and 'Sample
            %   Rate', with corresponding numeric values.
            %
            %   (3) Two of the column names must equal 'Trap1X' and 'Trap2X'. The
            %   corresponding columns will be treated as the data for beads A and B,
            %   respectively.
            %
            %   The outputs fs, KA, KB, CALA, and CALB are the sampling rate, stiffnesses
            %   for beads A and B, and volt-to-nanometer conversion factors for beads A
            %   and B, respectively.
            %
            %   The output tData is a table containing the data, with column names
            %   determined by MATLAB from colName1 through colNameM and with the columns
            %   for beads A and B renamed from 'Trap1X' and 'Trap2X' to 'BeadAPos' and
            %   'BeadBPos'. The data in these columns are in units of volts.
            %
            %   The output header is a 1x2 cell array. header{1} is a Nx1 cell array
            %   containing the keys, as character vectors. header{2} is a Nx1 cell array
            %   containing the values, also as character vectors.
            
            % Try to open the file.
            fileID = fopen(file);
            if fileID == -1
                error(['Unable to open ' file ' for reading.'])
            end
            
            try
                % Count the number of lines which contain non-numeric data.
                i = 0;
                while ~feof(fileID)
                    % While there are lines left to read, read the next line and determine the
                    % tab-delimited text elements.
                    C = regexp(fgetl(fileID), '[^\t]+', 'match');

                    % If all of the tab-delimited elements are numeric, exit the loop.
                    if all(~isnan(str2double(C))), break; end

                    % Otherwise, increment i.
                    i = i + 1;
                end

                % Read the header.
                frewind(fileID); % Reset the current line, so that textscan reads from line 1.
                header = textscan(fileID, '%s%s', i-1, 'Delimiter', '\t'); % Read tab-delimited key value pairs up until the column names.
                fclose(fileID);
            catch ME
                fclose(fileID);
                rethrow(ME);
            end
            
            % Read fs, KA, KB, CALA, and CALB from the header.
            KA = str2double(header{2}{strcmp(header{1}, 'K1')});
            KB = str2double(header{2}{strcmp(header{1}, 'K3')});
            CALA = str2double(header{2}{strcmp(header{1}, 'CAL1')});
            CALB = str2double(header{2}{strcmp(header{1}, 'CAL3')});
            fs = str2double(header{2}{strcmp(header{1}, 'Sample Rate')});
            
            % Read the data into a table.
            tData = readtable(file, 'HeaderLines', i-1);
            
            % Rename the columns that contain the data.
            tData.Properties.VariableNames{'Trap1X'} = 'BeadAPos';
            tData.Properties.VariableNames{'Trap2X'} = 'BeadBPos';
        end
        
        function installed = checkInstall(toolbox)
            % CHECKINSTALL Checks if the specified toolbox is installed.
            %   installed = CHECKINSTALL(toolbox) returns true if the specified toolbox is
            %   installed for the current version of MATLAB.
            %
            %   For example, CHECKINSTALL('Signal Processing Toolbox') returns true
            %   (logical 1) if the Signal Processing Toolbox is available and false
            %   (logical 0) otherwise.
            %
            % See also ver.
            
            if isdeployed
                installed = true;
            else
                v = ver;
                installed = any(strcmp({v.Name}, toolbox));
            end
        end
        
        function cov = calculateCov(A, B, averagingWindow, smoothingWindow, tlbxAvailable)
            % CALCULATECOV Calculates the covariance between two series of data.
            %   cov = CALCULATECOV(A, B, averagingWindow) calculates the
            %   sliding window covariance between A and B with window size determined by
            %   averagingWindow.
            %
            %   cov = CALCULATECOV(A, B, averagingWindow, smoothingWindow) passes the
            %   covariance through a second order Savitzky-Golay filter with window size
            %   smoothingWindow, provided that the Signal Processing Toolbox is installed.
            %   The function will check if this toolbox is installed.
            %
            %   cov = CALCULATECOV(A, B, averagingWindow, smoothingWindow, tlbxAvailable)
            %   does not check if the Signal Processing Toolbox is installed and instead
            %   trusts the value of tlbxAvailable (true or false). An error will occur if
            %   tlbxAvailable is set to true but the toolbox is not actually installed.
            
            if nargin < 5
                tlbxAvailable = SPASM.checkInstall('Signal Processing Toolbox');
            end
            
            % Calculate the covariance.
            A_filt = movmean(A, averagingWindow);
            B_filt = movmean(B, averagingWindow);
            AB_filt = movmean(A.*B, averagingWindow);
            cov = AB_filt - A_filt.*B_filt;
            % REQUIRES SIGNAL PROCESSING TOOLBOX
            if tlbxAvailable && numel(cov) > 2 && nargin > 3
                cov = sgolayfilt(cov, 2, smoothingWindow);
            end
        end
        
        function [p1, p2, m] = calculatePeaksMin(edges, values)
            % CALCULATEPEAKSMIN Calculates the two peaks of a bimodal histogram, as well
            % as the minimum value between the two peaks.
            %   [p1, p2, m] = CALCULATEPEAKSMIN(edges, values) returns the lower peak p1,
            %   the upper peak p2, and the minimum m of the histogram defined by inputs
            %   edges and values.
            %
            %   See also histcounts.
            
            % X and Y coordinates of tops of histogram bins.
            x = edges(1:end-1) + diff(edges) / 2;
            y = values;
            
            % X and Y coordinates of interpolated spline function.
            xx = linspace(x(1), x(end), 10*length(x));
            yy = spline(x, y, xx);
            
            [~, P] = islocalmax(yy); % [~, score indicating likelihood of each point being a local peak]
            [~, i] = maxk(P, 2); % [~, indices of top two scores]
            i = sort(i);
            [~, ii] = min(yy(i(1):i(2))); % [~, indices of minimum value in between top two scores]
            
            m = xx(ii+i(1)-1); % Minimum.
            p1 = xx(i(1)); % Peak 1.
            p2 = xx(i(2)); % Peak 2.
        end
        
        function prelimEvents = findPrelimEvents(cov, p1, p2)
            % FINDPRELIMEVENTS Finds events using the peak to peak method.
            %   prelimEvents = FINDPRELIMEVENTS(cov, p1, p2) uses the covariance cov and
            %   the peak locations p1 and p2 to find the indices within cov which mark the
            %   beginnings and ends of events as determined by the peak to peak method.
            %   These indices are returned in a 2-column array: the first column gives the
            %   event starts (the last index before each event), the second column gives
            %   the event ends (the last index of each event), and the rows are ordered
            %   chronologically.
            %
            % In the peak to peak method, event starts are found when the covariance drops
            % below the upper peak, provided it then drops below the lower peak. Event
            % ends are found when the covariance then rises above the upper peak again.
            
            startp2 = [diff(cov < p2); 0] == 1; % Find all points where the covariance drops below peak 2...
            istartp2 = find(startp2);

            stopp2 = [diff(cov < p2); 0] == -1; % ...and all points where the covariance rises above peak 2...
            stopp2(1:istartp2(1)) = 0; % ...excluding any such points without a preceding drop as found above...
            istopp2 = find(stopp2);

            startp1 = [diff(cov < p1); 0] == 1; % ...and find all points where the covariance drops below peak 1.
            istartp1 = find(startp1);

            % Consider events to occur when the covariance first drops below peak
            % 2, then drops below peak 1, and finally rises above peak 2.
            d = diff(sign(istartp2 - istartp1'));
            [ii, ~] = ind2sub(size(d), find(d > 0));
            ii = unique(ii);
            istart = istartp2(ii);
            istop = istopp2(ii);
            prelimEvents = [istart istop];
        end
        
        function prelimEvents = findPrelimEventsByMin(cov, m)
            % FINDPRELIMEVENTSBYMIN Finds events using the single threshold method.
            %   prelimEvents = FINDPRELIMEVENTSBYMIN(cov, m) uses the covariance cov and
            %   the minimum location m to find the indices within cov which mark the
            %   beginnings and ends of events as determined by the single threshold
            %   method. These indices are returned in a 2-column array: the first column
            %   gives the event starts (the last index before each event), the second
            %   column gives the event ends (the last index of each event), and the rows
            %   are ordered chronologically.
            %
            % In the single threshold method, event starts are found when the covariance
            % drops below the minimum, and event ends are found when the covariance rises
            % above the minimum.
            
            index = cov < m; % Select areas where the covariance is less than the threshold...
            start = [diff(index); 0] == 1; % ...find the beginnings of these areas...
            stop = [diff(index); 0] == -1; % ...and the ends of these areas...
            stop(find(stop, 1)) = find(start, 1) < find(stop, 1); % ...remove any areas which start before file starts...
            start(find(start, 1, 'last')) = find(start, 1, 'last') < find(stop, 1, 'last'); % ...as well as any areas which end after file ends...
            istart = find(start); % ...and store the indices where these areas begin...
            istop = find(stop); % ...and where they end.
            prelimEvents = [istart istop];
        end
        
        function realEvents = findRealEvents(key)
            % FINDREALEVENTS Finds events according to a key.
            %   realEvents = FINDREALEVENTS(key) finds the indices within column vector
            %   key which mark the beginnings and ends of events. These indices are
            %   returned in a 2-column array: the first column gives the event starts (the
            %   last index before each event), the second column gives the event ends (the
            %   last index of each event), and the rows are ordered chronologically. Any
            %   non-zero element within key is considered part of an event.
            %
            % For example, with input key = [0; 0; 1; 1; 1; 0; 1; 1; 0], the output would
            % be realEvents = [3, 6; 7, 9].
            
            key(key ~= 0) = 1; % Any non-zero element is considered part of a binding event.
            start = [diff(key); 0] == 1; % Logical array marking event starts.
            stop = [diff(key); 0] == -1; % Logical array marking event ends.
            stop(find(stop, 1)) = find(start, 1) < find(stop, 1); % Remove events starting before start of data.
            start(find(start, 1, 'last')) = find(start, 1, 'last') < find(stop, 1, 'last'); % Remove events ending after end of data.
            realEvents = [find(start), find(stop)];
        end
        
        function [FP_FN, startError, endError] = scoreDetectedEvents(detected, real)
            % SCOREDETECTEDEVENTS Compares detected events to real events.
            %   FP_FN = SCOREDETECTEDEVENTS(detected, real) compares detected
            %   events in detected to real events in real, where detected and real are
            %   2-column arrays with event starts in the first column and event ends in
            %   the second column.
            %   
            %   Each detected event is mapped to the closest overlapping real event. Each
            %   real event is mapped to the closest overlapping detected event. Events
            %   which are mapped to each other are considered correctly identified events.
            %   All other detected events are considered false positives. All other real
            %   events are considered false negatives.
            %
            %   The struct FP_FN contains the following fields:
            %     FP - the number of false positive events
            %     FN - the number of false negative events
            %     detectedIdx - the indices of the detected events which were NOT
            %       counted as false positive events
            %     realIdx - the indices of the real events which were NOT counted as
            %       false negative events
            %
            %   [FP_FN, startError, endError] = SCOREDETECTEDEVENTS(detected, real) also
            %   compares the start and end of each detected event to the start and end of
            %   each real event, respectively, to find the error.
            %
            %   Each detected event's start is mapped to the closest start of all
            %   overlapping real events. Each real event's start is mapped to the closest
            %   start of all overlapping detected events. Starts which are mapped to each
            %   other are compared to find the error in the event starts.
            %
            %   The error in the event ends is found similarly.
            %
            %   The struct startError contains the following fields:
            %     error - the error in the event starts
            %     detectedIdx - the indices of the detected events whose starts were
            %       considered in the calculation of error
            %     realIdx - the indices of the real events whose starts were considered in
            %       the calculation of error
            %
            %   The struct endError contains similar fields.
            
            DtoR = nan(size(detected,1),1); % Detected event i maps to real event DtoR(i).
            RtoD = nan(size(real,1),1); % Real event i maps to detected event RtoD(i).
            sDtoR = DtoR; % The start of detected event i maps to the start of real event sDtoR(i).
            eDtoR = DtoR; % The end of detected event i maps to the end of real event sDtoR(i).
            sRtoD = RtoD; % The start of real event i maps to the start of detected event sDtoR(i).
            eRtoD = RtoD; % The end of real event i maps to the end of detected event sDtoR(i).
            
            for i = 1:size(detected,1)
                idx = find(real(:,2) > detected(i,1) & real(:,1) < detected(i,2));
                if ~isempty(idx)
                    dStart = abs(detected(i,1)-real(idx,1));
                    dEnd = abs(detected(i,2)-real(idx,2));
                    [~, minIdx] = min((dStart + dEnd)/2);
                    DtoR(i) = idx(minIdx);
                    [~, sMinIdx] = min(dStart);
                    sDtoR(i) = idx(sMinIdx);
                    [~, eMinIdx] = min(dEnd);
                    eDtoR(i) = idx(eMinIdx);
                end
            end
            
            for i = 1:size(real,1)
                idx = find(detected(:,2) > real(i,1) & detected(:,1) < real(i,2));
                if ~isempty(idx)
                    dStart = abs(real(i,1)-detected(idx,1));
                    dEnd = abs(real(i,2)-detected(idx,2));
                    [~, minIdx] = min((dStart + dEnd)/2);
                    RtoD(i) = idx(minIdx);
                    [~, sMinIdx] = min(dStart);
                    sRtoD(i) = idx(sMinIdx);
                    [~, eMinIdx] = min(dEnd);
                    eRtoD(i) = idx(eMinIdx);
                end
            end
            
            % Indices (within real and detected) of correctly identified events:
            A = RtoD == (1:numel(DtoR));
            B = DtoR == (1:numel(RtoD));
            [FP_FN.realIdx, FP_FN.detectedIdx] = find(A & B');
            % Real event realIdx(i) and detected event detectedIdx(i) were mapped to each
            % other.
            
            % Events not mapped to each other are either false positive or false negative
            % events.
            FP_FN.FN = size(real,1) - numel(FP_FN.realIdx);
            FP_FN.FP = size(detected,1) - numel(FP_FN.detectedIdx);
            
            % Only event starts which mapped to each other contribute to the error.
            A = sRtoD == (1:numel(sDtoR));
            B = sDtoR == (1:numel(sRtoD));
            [startError.realIdx, startError.detectedIdx] = find(A & B');
            
            % Only event ends which mapped to each other contribute to the error.
            A = eRtoD == (1:numel(eDtoR));
            B = eDtoR == (1:numel(eRtoD));
            [endError.realIdx, endError.detectedIdx] = find(A & B');
            
            startError.error = detected(startError.detectedIdx,1) - real(startError.realIdx,1);
            endError.error = detected(endError.detectedIdx,2) - real(endError.realIdx,2);
        end
        
        function [numPtsBefore, numPtsAfter] = findChangepointWindows(events, N)
            % FINDCHANGEPOINTWINDOWS Given a collection of events, finds windows
            % surrounding each event for the changepoint algorithm to consider.
            %
            %   [numPtsBefore, numPtsAfter] = FINDCHANGEPOINTWINDOWS(events, N) returns,
            %   for each event in events, the number of points before the start of the
            %   event and the number of points after the end of the event. events should
            %   be a 2-column array: the first column should give the event starts (the
            %   last index before each event), the second column should give the event
            %   ends (the last index of each event), and the rows should be ordered
            %   chronologically. N is the total number of data points, which is needed to
            %   determine the window for the final event.
            
            events = [0 0; events; N, N];
            sep = events(2:end,1) - events(1:end-1,2);
            dur = events(2:end-1,2) - events(2:end-1,1);
            numPtsBefore = floor(0.49*min(sep(1:end-1), dur));
            numPtsAfter = floor(0.49*min(sep(2:end), dur));
        end
        
        function [T1, T2, Lr, idx] = changepoint(data, varargin)
            % CHANGEPOINT Finds two changepoints in a collection of points.
            %   [T1, T2] = CHANGEPOINT(data) returns the indices in data, T1 and T2, which
            %   are most likely the indices of the true changepoints. It is assumed that
            %   two true changepoints exist such that all points between these true
            %   changepoints are drawn from one normal distribution and all other points
            %   are drawn from another normal distribution. The variance between the two
            %   changepoints will be smaller than the variance before the first
            %   changepoint and the variance after the second changepoint. T1 is the index
            %   of the last point of the first section of data (i.e. the last point before
            %   the event). T2 is the index of the last point of the second section of
            %   data (i.e. the last point of the event).
            %
            %   [T1, T2, Lr, idx] = CHANGEPOINT(data) also returns the likelihood Lr.
            %   Lr(i,j) gives the relative likelihood that points idx(i) and idx(j) are
            %   the true changepoints. Note that for events longer than 1500 points, this
            %   method does not calculate the likelihood for every single pair of points.
            %   Rather, the likelihood is calculated for pairs of points drawn from 1500
            %   evenly spaced points throughout the data, providing initial estimates of
            %   T1 and T2. Then, for each estimate, this method zooms in on the data
            %   surrounding that estimate and runs the changepoint algorithm again. idx
            %   gives the indices of the 1500 evenly spaced points.
            %
            %   [___] = CHANGEPOINT(data, 'likelihood') only calculates and returns Lr and
            %   idx. T1 and T2 will be empty.
            
            onlyLr = any(strcmp(varargin, 'likelihood'));
            
            data = data(:);
            if numel(data) <= 1500
                % If the event is short, find two changepoints within the data and return.
                [T1, T2, Lr] = SPASM.findTwoChangepoints(data, onlyLr);
                idx = 1:numel(data);
            else
                % If the event is too long, consider a subset of the data pulled from 1500
                % evenly spaced points.
                idx = floor(linspace(1, numel(data), 1500));
                subset = data(idx);
                [t1, t2, Lr] = SPASM.findTwoChangepoints(subset, onlyLr); % [T1's estimate, T2's estimate, likelihood]
                
                if onlyLr
                    % If only the likelihood was requested, return.
                    T1 = t1;
                    T2 = t2;
                    return
                end
                
                % Zoom in around T1's estimate to find T1.
                idx1 = max(idx(t1)-1000,1):min(idx(t1)+1000,max(idx(t2)-10,idx(t1)+1));
                subset1 = data(idx1);
                T1 = SPASM.findOneChangepoint(subset1);
                T1 = idx1(T1);
                
                % Zoom in around T2's estimate to find T2.
                idx2 = max(idx(t2)-1000,min(idx(t1)+10,idx(t2)-1)):min(idx(t2)+1000,numel(data));
                subset2 = data(idx2);
                T2 = SPASM.findOneChangepoint(subset2);
                T2 = idx2(T2);
            end
        end
        
        function [T1, T2, Lr] = findTwoChangepoints(data, onlyLr)
            % FINDTWOCHANGEPOINTS Helper function for changepoint() which finds two
            % changepoints in a set of data.
            %   [T1, T2] = FINDTWOCHANGEPOINTS(data) returns the indices in data, T1 and
            %   T2, which are most likely the indices of the true changepoints. It is
            %   assumed that two true changepoints exist such that all points between
            %   these true changepoints are drawn from one normal distribution and all
            %   other points are drawn from another normal distribution. The variance
            %   between the two changepoints will be smaller than the variance before the
            %   first changepoint and the variance after the second changepoint. T1 is the
            %   index of the last point of the first section of data (i.e. the last point
            %   before the event). T2 is the index of the last point of the second section
            %   of data (i.e. the last point of the event).
            %
            %   [T1, T2, Lr] = FINDTWOCHANGEPOINTS(data) also returns the likelihood Lr.
            %   Lr(i,j) gives the relative likelihood that points i and j are the true
            %   changepoints. This method calculates the likelihood for every single pair
            %   of points, regardless of the length of data.
            %
            %   [___] = FINDTWOCHANGEPOINTS(data, onlyLr) does not bother calculating T1
            %   and T2 and instead returns empty arrays for both as long as onlyLr is true
            %   (logical 1). onlyLr is false by default.
            %
            % If data is very long, this method will use a lot of memory and may crash
            % MATLAB.
            %
            % See also changepoint.
            %
            % Created with help from:
            % Killick R, Haynes K and Eckley IA (2016). changepoint: An R
            % package for changepoint analysis. R package version 2.2.2,
            % https://CRAN.R-project.org/package=changepoint.
            
            if nargin == 1
                onlyLr = false;
            end
            data = data(:);
            N = length(data);
            if N < 5
                T1 = 1;
                T2 = N;
                Lr = zeros(N);
                return
            end
            tC = (1:N)'; % Vector of candidate points.
            y = cumsum(data);
            y2 = cumsum(data.^2);
            s20 = (y2(N) - (y(N)^2/N))/N; % Variance of entire data set.
            s213s = ((y2(N)-y2+y2') - ((y(N)-y+y').^2./(N-tC+tC')))./(N-tC+tC'); % For each pair of candidate points t1 and t2, variance of all points between t1 and t2.
            s213s(s213s < 0) = 0; % (For those numbers which are essentially zero yet negative for some reason.)
            s22s = ((y2-y2') - ((y-y').^2./(tC-tC')))./(tC-tC'); % For each pair of candidate points t1 and t2, variance of all points between t1 and t2.
            s22s(s22s < 0) = 0;
            
            % Calculate the likelihood that each pair of candidate points is the pair of
            % true changepoints. Lr(t1,t2) gives the likelihood that t1 and t2 are the
            % true changepoints.
            Lr = (N/2)*log(s20) - ((tC'-tC)/2).*log(s22s) - ((N-tC'+tC)/2).*log(s213s);
            Lr(1,:) = 0; % First changepoint should be at least 2...
            Lr(end-3:end, :) = 0; % ...and no greater than N-4.
            Lr(:,1:3) = 0; % Second changepoint should be at least 4...
            Lr(:,end-1:end) = 0; % ...and no greater than N-2.
            Lr = triu(Lr, 2); % Second changepoint should be at least 2 greater than the first.
            Lr(isinf(Lr)) = NaN; % Ignore infinite values.
            
            % If only Lr is requested, no need to go any further.
            if onlyLr
                T1 = [];
                T2 = [];
                return
            end
            
            % Pick best pair among the pairs for which variance drops during the event.
            validIdx = s213s > s22s;
            [~, maxIdx] = max(Lr(validIdx));
            [row, col] = find(validIdx, maxIdx);
            row = row(end);
            col = col(end);
            T1 = min([row, col]);
            T2 = max([row, col]);
        end
        
        function [T, Lr] = findOneChangepoint(data)
            % FINDONECHANGEPOINT finds a single changepoint in data.
            %   T = FINDONECHANGEPOINT(data) returns the index in data, T, which is most
            %   likely the index of the true changepoint. It is assumed that one true
            %   changepoint exists such that all points before and including this true
            %   changepoint are drawn from one normal distribution and all other points
            %   are drawn from another normal distribution.
            %
            %   [T, Lr] = FINDONECHANGEPOINT(data) also returns the likelihood Lr. Lr(i)
            %   gives the relative likelihood that point i is the true changepoint. This
            %   method calculates the likelihood for every single point, regardless of the
            %   length of data.
            %
            % Created with help from:
            % Killick R, Haynes K and Eckley IA (2016). changepoint: An R
            % package for changepoint analysis. R package version 2.2.2,
            % https://CRAN.R-project.org/package=changepoint.
            
            data = data(:);
            N = length(data);
            if N < 3
                T = 1;
                Lr = zeros(N);
                return
            end
            tC = (1:N)'; % Vector of candidate points.
            y = cumsum(data);
            y2 = cumsum(data.^2);
            s20 = (y2(N) - (y(N)^2/N))/N; % Variance of entire data set.
            s21s = (y2 - (y.^2./tC))./tC; % For each candidate point t1, variance of all points before t1.
            s21s(s21s < 0) = 0; % (For those numbers which are essentially zero yet negative for some reason.)
            s22s = ((y2(N)-y2) - ((y(N)-y).^2./(N-tC)))./(N-tC); % For each candidate point t1, variance of all points after t1.
            s22s(s22s < 0) = 0;

            % Calculate the likelihood that each pair of candidate points is the pair of
            % true changepoints. Lr(t1,t2) gives the likelihood that t1 and t2 are the
            % true changepoints.
            Lr = (N/2)*log(s20) - (tC/2).*log(s21s) - (((N-tC)/2).*log(s22s));
            Lr(1) = 0; % Changepoint should be at least 2...
            Lr(end-1:end) = 0; % ...and no greater than N-2.
            Lr(isinf(Lr)) = NaN; % Ignore infinite values.

            [~, T] = max(Lr); % [~, index of maximum likelihood]
        end
    end
    
    % Private general methods
    methods (Access = private)
        
        function initialize(app)
            % Initializes all properties.
            
            % Colors:
            app.misc.colors.beadA = [0, 0.447, 0.741];
            app.misc.colors.beadB = [0.85, 0.325, 0.098];
            app.misc.colors.mainPatches = [1, 0, 0];
            app.misc.colors.mainPatchesSelected = [0.9290, 0.6940, 0.1250];
            app.misc.colors.simPatches = [0, 0, 1];
            app.misc.colors.deselect = [1, 0, 0];
            app.misc.colors.deselectSelected = [0, 0, 1];
            app.misc.colors.cov = [0, 0.447, 0.741];
            app.misc.colors.covPatches = [1, 0, 0];
            app.misc.colors.covHist = [0, 0.447, 0.741];
            app.misc.colors.peaksMin = [1, 0, 0];
            app.misc.colors.windowsPt = [0, 0.447, 0.741];
            app.misc.colors.misc = [0, 0.447, 0.741];
            app.misc.colors.miscFit = [1, 0, 0];
            app.misc.colors.ensStart = [0, 0.447, 0.741];
            app.misc.colors.ensEnd = [0, 0.447, 0.741];
            app.misc.colors.ensStartFit = [1, 0, 0];
            app.misc.colors.ensEndFit = [1, 0, 0];
            app.misc.colors.individual = [0, 0.447, 0.741];
            app.misc.colors.defChngpt = [0, 0, 0];
            app.misc.colors.chngpt = [1, 0, 0];
            app.misc.colors.simLines = [0.9290, 0.6940, 0.1250];
            app.misc.colors.removeLabel = [1, 0, 0];
            app.misc.colors.includeLabel = [0.1, 0.7, 0.1];
            app.misc.colors.currentTask = [1, 0, 0];
            app.misc.colors.activeAxes = [1, 0, 0];
            
            % Figure:
            screen = get(0, 'ScreenSize');
            width = 0.83 * screen(3);
            height = 0.83 * screen(4);
            app.misc.fontSize = min(height / 101.25, width / 180);
            app.misc.fig = figure('Visible', 'off',...
                'Position', [0, 0, width, height],...
                'MenuBar', 'none',...
                'Toolbar', 'none',...
                'Resize', 'on',...
                'Name', 'SPASM',...
                'NumberTitle', 'off',...
                'WindowButtonMotionFcn', @app.defaultMotion,...
                'WindowButtonDownFcn', @app.defaultClick,...
                'WindowKeyPressFcn', @app.nullKeyPressCallback,...
                'CloseRequestFcn', @app.closeRequest);
            app.misc.storeClickCallback = app.misc.fig.WindowButtonDownFcn;
            app.misc.storeMotionCallback = app.misc.fig.WindowButtonMotionFcn;
            app.misc.cm = uicontextmenu(app.misc.fig);
            
            % Tabs:
            app.misc.tabgroup = uitabgroup(app.misc.fig);
            app.misc.tab1 = uitab(app.misc.tabgroup,...
                'Title', 'Analyze Single File');
            app.misc.tab2 = uitab(app.misc.tabgroup,...
                'Title', 'Combine Multiple Files');
            
            % Load panel of tab 1:
            app.load.panel = uipanel(app.misc.tab1,...
                'Title', 'Load',...
                'Position', [0.005, 0.52, 0.99, 0.48],...
                'FontSize', app.misc.fontSize);
            app.load.inputBtn = uicontrol(app.load.panel,...
                'Style', 'pushbutton',...
                'String', 'Input Raw Data',...
                'Units', 'normalized',...
                'Position', [0.005, 0.902, 0.085, 0.083],...
                'FontSize', app.misc.fontSize,...
                'Callback', @app.chooseTxtCallback);
            app.load.txtInputFcn = @app.defaultTxtInput;
            app.load.txtOutputFcn = @app.defaultTxtOutput;
            app.load.path = [];
            app.load.name = [];
            app.load.ext = [];
            app.load.header = [];
            app.data.origData = [];
            app.data.origBeads = [];
            app.data.time = [];
            app.load.CALA = [];
            app.load.CALATxt = uicontrol(app.load.panel,...
                'Units', 'normalized',...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load.CALB = [];
            app.load.CALBTxt = uicontrol(app.load.panel,...
                'Units', 'normalized',...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load.KA = [];
            app.load.KATxt = uicontrol(app.load.panel,...
                'Units', 'normalized',...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load.KB = [];
            app.load.KBTxt = uicontrol(app.load.panel,...
                'Units', 'normalized',...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.data.fs = [];
            app.load.fsTxt = uicontrol(app.load.panel,...
                'Units', 'normalized',...
                'Style', 'text',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load.filenameTxt = uicontrol(app.load.panel,...
                'Units', 'normalized',...
                'Style', 'text',...
                'String', '',...
                'HorizontalAlignment', 'left',...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.misc.simInput = false;
            app.misc.realEvents = [];
            app.misc.simLines = gobjects(0);
            app.misc.simPatches = gobjects(0);
            
            app.load.currentTask = uicontrol(app.load.panel,...
                'Style', 'text',...
                'String', '',...
                'ForegroundColor', app.misc.colors.currentTask,...
                'Units', 'normalized',...
                'Position', [0.99 0.93 0.01 0.055],...
                'FontSize', 1.3*app.misc.fontSize,...
                'HorizontalAlignment', 'right');
            app.load.taskMaxWidth = 1;
            app.load.skipBtn = uicontrol(app.load.panel,...
                'Style', 'pushbutton',...
                'String', 'Skip',...
                'Units', 'normalized',...
                'FontSize', app.misc.fontSize,...
                'Callback', @app.skipBtnCallback,...
                'Visible', 'off');
            app.misc.skip = false;
            
            % Main axes of load panel of tab 1:
            app.load.main.mainAxes = axes(app.load.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.005, 0.352, 0.965, 0.532],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'mainAxes');
            xlabel(app.load.main.mainAxes, 'time (s)');
            ylabel(app.load.main.mainAxes, 'position (nm)');
            app.load.main.ylim = [];
            app.load.main.bead1Plot = plot(app.load.main.mainAxes, 0, 0,...
                'Color', app.misc.colors.beadA,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.load.main.bead2Plot = plot(app.load.main.mainAxes, 0, 0,...
                'Color', app.misc.colors.beadB,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.load.main.patches = gobjects(0);
            app.load.main.label = annotation(app.load.panel, 'textbox',...
                'FitBoxToText', 'on',...
                'BackgroundColor', 'white',...
                'Visible', 'off',...
                'Units', 'pixels');
            app.load.main.prevIn = [];
            
            % Covariance axes of load panel of tab 1:
            app.load.cov.covAxes = axes(app.load.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.005, 0.092, 0.965, 0.25],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'covAxes');
            xlabel(app.load.cov.covAxes, 'time (s)');
            ylabel(app.load.cov.covAxes, 'covariance');            
            app.load.cov.covPlot = plot(app.load.cov.covAxes, 0, 0,...
                'Color', app.misc.colors.cov,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.load.cov.peak1Plot = plot(app.load.cov.covAxes, 0, 0,...
                'Color', app.misc.colors.peaksMin,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.load.cov.peak2Plot = plot(app.load.cov.covAxes, 0, 0,...
                'Color', app.misc.colors.peaksMin,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.load.cov.minPlot = plot(app.load.cov.covAxes, 0, 0,...
                'Color', app.misc.colors.peaksMin,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.load.cov.patches = gobjects(0);
            app.load.cov.cov = [];
            
            % Analyze panel of tab 1:
            app.analyze.panel = uipanel(app.misc.tab1,...
                'Title', 'Analyze',...
                'Position', [0.005, 0.005, 0.6, 0.513],...
                'FontSize', app.misc.fontSize);
            
            % Covariance histogram axes of analyze panel of tab 1:
            app.analyze.covHist.covHistAxes = axes(app.analyze.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.01, 0.581, 0.284, 0.397],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'covHistAxes');
            xlabel(app.analyze.covHist.covHistAxes, 'covariance');
            ylabel(app.analyze.covHist.covHistAxes, 'count');
            app.analyze.covHist.covHistPlot = histogram(app.analyze.covHist.covHistAxes, [],...
                'FaceColor', app.misc.colors.covHist,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.covHist.peak1Plot = plot(app.analyze.covHist.covHistAxes, 0, 0,...
                'Color', app.misc.colors.peaksMin,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.covHist.peak2Plot = plot(app.analyze.covHist.covHistAxes, 0, 0,...
                'Color', app.misc.colors.peaksMin,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.covHist.minPlot = plot(app.analyze.covHist.covHistAxes, 0, 0,...
                'Color', app.misc.colors.peaksMin,...
                'Visible', 'off',...
                'HitTest', 'off');
            
            % Miscellaneous axes of analyze panel of tab 1:
            app.analyze.misc.miscAxes = axes(app.analyze.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.34, 0.581, 0.284, 0.397],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'miscAxes');
            app.analyze.misc.miscPlot = plot(app.analyze.misc.miscAxes, 0, 0,...
                'Color', app.misc.colors.misc,...
                'MarkerSize', 3,...
                'MarkerFaceColor', app.misc.colors.misc,...
                'MarkerEdgeColor', app.misc.colors.misc,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.misc.realPlot = plot(app.analyze.misc.miscAxes, 0, 0,...
                'Color', app.misc.colors.simLines,...
                'MarkerSize', 3,...
                'MarkerFaceColor', app.misc.colors.simLines,...
                'MarkerEdgeColor', app.misc.colors.simLines,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.misc.miscFitPlot = plot(app.analyze.misc.miscAxes, 0, 0,...
                'Color', app.misc.colors.miscFit,...
                'LineWidth', 1.5,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.misc.miscLabel = text(app.analyze.misc.miscAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.misc.data = 'Event Duration'; % Show event durations by default.
            app.analyze.misc.kDur = [];
            app.analyze.misc.meanStep1 = [];
            app.analyze.misc.stdStep1 = [];
            app.analyze.misc.meanStep2 = [];
            app.analyze.misc.stdStep2 = [];
            app.analyze.misc.meanTotalStep = [];
            app.analyze.misc.stdTotalStep = [];
            
            % Ensemble average axes of analyze panel of tab 1:
            app.analyze.ensemble.ensembleAxes = axes(app.analyze.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.67, 0.581, 0.284, 0.397],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'ensembleAxes');
            xlabel(app.analyze.ensemble.ensembleAxes, 'time (s)');
            app.analyze.ensemble.globalFit = true;
            app.analyze.ensemble.start = [];
            app.analyze.ensemble.end = [];
            app.analyze.ensemble.startFit = [];
            app.analyze.ensemble.endFit = [];
            app.analyze.ensemble.p = [];
            app.analyze.ensemble.startPlot = plot(app.analyze.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensStart,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.ensemble.endPlot = plot(app.analyze.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensEnd,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.ensemble.startFitPlot = plot(app.analyze.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensStartFit,...
                'LineWidth', 2,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.ensemble.endFitPlot = plot(app.analyze.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensEndFit,...
                'LineWidth', 2,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.ensemble.startLabel = text(app.analyze.ensemble.ensembleAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.ensemble.endLabel = text(app.analyze.ensemble.ensembleAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.ensemble.stepLabel = text(app.analyze.ensemble.ensembleAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.ensemble.pts_bef = 0.75; % seconds
            app.analyze.ensemble.pts_after = 0.75; % seconds
            app.analyze.ensemble.exten_pos = 0.01; % seconds
            app.analyze.ensemble.num_pts_to_skip = 10; % number of points from transition to skip when fitting
            
            % Likelihood axes of analyze panel of tab 1:
            app.analyze.individual.pickChngptBtn = uicontrol(app.analyze.panel,...
                'Style', 'pushbutton',...
                'String', '<html><p align="center">Pick<br/>Manually</p></html>',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'center',...
                'Position', [0.007, 0.324, 0.07, 0.165],...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'Callback', @app.pickChngptCallback);
            app.analyze.individual.resetBtn = uicontrol(app.analyze.panel,...
                'Style', 'pushbutton',...
                'String', '<html><p align="center">Reset to<br/>Default</p></html>',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'center',...
                'Position', [0.007, 0.154, 0.07, 0.165],...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'Callback', @app.defChngptCallback);
            app.analyze.individual.snap = uicontrol(app.analyze.panel,...
                'Style', 'checkbox',...
                'String', '<html><p align="center">Snap to<br/>Change Point</p></html>',...
                'Units', 'normalized',...
                'Position', [0.007, 0.029, 0.07, 0.12],...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off');
            app.analyze.individual.likelihoodAxes = axes(app.analyze.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.084, 0.092, 0.19, 0.397],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YTick', [],...
                'YLabel', [],...
                'ZTick', [],...
                'ZLabel', [],...
                'UserData', 'likelihoodAxes');
            xlabel(app.analyze.individual.likelihoodAxes, 'time (s)');
            app.analyze.individual.surfPlot = surf(app.analyze.individual.likelihoodAxes, [], [], [],...
                'Visible', 'off',...
                'EdgeColor', 'none',...
                'HitTest', 'off');
            app.analyze.individual.defXonPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defXoffPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defYonPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defYoffPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.XonPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.XoffPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.YonPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.YoffPlot = plot3(app.analyze.individual.likelihoodAxes, 0, 0, 0, '--',...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.data.defaultEvents = [];
            app.data.windowStarts = [];
            app.data.windowStops = [];
            app.data.defaultWindowStarts = [];
            app.data.defaultWindowStops = [];
            
            % Individual event axes of analyze panel of tab 1:
            app.analyze.individual.eventAxes = axes(app.analyze.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.281, 0.092, 0.434, 0.397],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'eventAxes');
            xlabel(app.analyze.individual.eventAxes, 'time (s)');
            app.analyze.individual.eventPlot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.individual,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defTonPlot = plot(app.analyze.individual.eventAxes, 0, 0, '--',...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defToffPlot = plot(app.analyze.individual.eventAxes, 0, 0, '--',...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.TonPlot = plot(app.analyze.individual.eventAxes, 0, 0, '--',...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.ToffPlot = plot(app.analyze.individual.eventAxes, 0, 0, '--',...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.first = true;
            app.analyze.individual.T1 = [];
            app.analyze.individual.defMean1Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defUpLim1Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defLowLim1Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defMean2Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defUpLim2Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defLowLim2Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defMean3Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defUpLim3Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.defLowLim3Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.defChngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.mean1Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.upLim1Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.lowLim1Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.mean2Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.upLim2Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.lowLim2Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.mean3Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.upLim3Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze.individual.lowLim3Plot = plot(app.analyze.individual.eventAxes, 0, 0,...
                'Color', app.misc.colors.chngpt,...
                'Visible', 'off',...
                'HitTest', 'off');
            
            % Individual event panel of analyze panel of tab 1:
            app.analyze.individual.panel = uipanel(app.analyze.panel,...
                'Position', [0.755, 0.02, 0.235, 0.47],...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.totalPotentialTxt = uicontrol(app.analyze.individual.panel,...
                'Style', 'text',...
                'String', 'Number of Potential Events:',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'center',...
                'Position', [0, 0.85, 1, 0.1],...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.totalActiveTxt = uicontrol(app.analyze.individual.panel,...
                'Style', 'text',...
                'String', 'Number of Analyzed Events:',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'center',...
                'Position', [0, 0.73, 1, 0.1],...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.showRemovedBtn = uicontrol(app.analyze.individual.panel,...
                'Style', 'checkbox',...
                'String', 'Show Removed Events',...
                'Callback', @app.showRemovedCallback,...
                'Units', 'normalized',...
                'Position', [0.17, 0.62, 0.8, 0.1],...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.currentEvent = 1;
            app.analyze.individual.storeCurrentEvent = 1;
            app.analyze.individual.currentEventTxt = uicontrol(app.analyze.individual.panel,...
                'Style', 'text',...
                'String', 'Current Event:',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'left',...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.currentEventTxt.Position =...
                [0.17, 0.49, app.analyze.individual.currentEventTxt.Extent(3), 0.1];
            app.analyze.individual.currentEventEdit = uicontrol(app.analyze.individual.panel,...
                'Style', 'edit',...
                'String', '1',...
                'Units', 'normalized',...
                'Position', [app.analyze.individual.currentEventTxt.Extent(3)+0.18, 0.493, 0.15, 0.1],...
                'HorizontalAlignment', 'left',...
                'FontSize', app.misc.fontSize,...
                'Callback', @app.currentEventCallback);
            app.analyze.individual.ofNumTxt = uicontrol(app.analyze.individual.panel,...
                'Style', 'text',...
                'String', 'of ',...
                'Units', 'normalized',...
                'Position', [app.analyze.individual.currentEventTxt.Extent(3)+0.34, 0.49, 1, 0.1],...
                'HorizontalAlignment', 'left',...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.prevEventBtn = uicontrol(app.analyze.individual.panel,...
                'Style', 'pushbutton',...
                'String', 'Show Previous',...
                'Callback', @app.prevEventCallback,...
                'Units', 'normalized',...
                'Position', [0.05, 0.31, 0.425, 0.16],...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.nextEventBtn = uicontrol(app.analyze.individual.panel,...
                'Style', 'pushbutton',...
                'String', 'Show Next',...
                'Callback', @app.nextEventCallback,...
                'Units', 'normalized',...
                'Position', [0.525, 0.31, 0.425, 0.16],...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.removeBtn = uicontrol(app.analyze.individual.panel,...
                'Style', 'checkbox',...
                'String', 'Remove this Event',...
                'Callback', @app.removeEventCallback,...
                'Units', 'normalized',...
                'Position', [0.25, 0.05, 0.7, 0.1],...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.label = uicontrol(app.analyze.individual.panel,...
                'Style', 'text',...
                'String', '',...
                'Units', 'normalized',...
                'HorizontalAlignment', 'center',...
                'Position', [0.05, 0.18, 0.9, 0.1],...
                'FontSize', app.misc.fontSize);
            app.analyze.individual.removedByDur = [];
            app.analyze.individual.removedBySep = [];
            app.analyze.individual.removedByUser = [];
            app.analyze.individual.includedByUser = [];
            
            % Controls panel of tab 1:
            app.controls.panel = uipanel(app.misc.tab1,...
                'Title', 'l',...
                'Position', [0.605, 0.005, 0.39, 0.513],...
                'BorderType', 'none',...
                'FontSize', app.misc.fontSize);
            app.controls.panel.ForegroundColor = app.controls.panel.BackgroundColor;
            app.controls.generateCovBtn = uicontrol(app.controls.panel,...
                'Style', 'pushbutton',...
                'String', 'Generate Covariance',...
                'Callback', @app.generateCovCallback,...
                'Units', 'normalized',...
                'Position', [0.02, 0.79, 0.46, 0.1],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.findEventsBtn = uicontrol(app.controls.panel,...
                'Style', 'pushbutton',...
                'String', 'Find Events',...
                'Callback', @app.findEventsCallback,...
                'Units', 'normalized',...
                'Position', [0.02, 0.58, 0.46, 0.2],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.useChangepointBtn = uicontrol(app.controls.panel,...
                'Style', 'checkbox',...
                'String', 'Use the Change Point Algorithm',...
                'Value', 1,...
                'Callback', @app.useChangepointCallback,...
                'Enable', 'off',...
                'Units', 'normalized',...
                'Position', [0.1, 0.5, 0.35, 0.07],...
                'FontSize', app.misc.fontSize);
            app.controls.toggleEventsBtn = uicontrol(app.controls.panel,...
                'Style', 'checkbox',...
                'String', 'Show Events',...
                'Value', 1,...
                'Callback', @app.toggleEventsCallback,...
                'Enable', 'off',...
                'Units', 'normalized',...
                'Position', [0.07, 0.43, 0.22, 0.07],...
                'FontSize', app.misc.fontSize);
            app.controls.flipDataBtn = uicontrol(app.controls.panel,...
                'Style', 'checkbox',...
                'String', 'Flip Data',...
                'Callback', @app.flipDataCallback,...
                'Enable', 'off',...
                'Units', 'normalized',...
                'Position', [0.3, 0.43, 0.2, 0.07],...
                'FontSize', app.misc.fontSize);
            app.controls.eventDetectMethodBtnGrp = uibuttongroup(app.controls.panel,...
                'Position', [0.01, 0.32, 0.48, 0.1],...
                'SelectionChangedFcn', @app.detectionMethodCallback);
            app.controls.peakToPeakBtn = uicontrol(app.controls.eventDetectMethodBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Peak-to-Peak',...
                'Units', 'normalized',...
                'Position', [0.1, 0.05, 0.38, 0.9],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.thresholdBtn = uicontrol(app.controls.eventDetectMethodBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Single Threshold',...
                'Units', 'normalized',...
                'Position', [0.52, 0.05, 0.43, 0.9],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.peakDetectMethodBtnGrp = uibuttongroup(app.controls.panel,...
                'Position', [0.01, 0.11, 0.48, 0.2],...
                'SelectionChangedFcn', @app.manvautoCallback);
            app.controls.automaticBtn = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Auto',...
                'Units', 'normalized',...
                'Position', [0.05, 0.525, 0.25, 0.45],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.manualBtn = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Manual',...
                'Units', 'normalized',...
                'Position', [0.32, 0.525, 0.3, 0.45],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.setManualBtn = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'pushbutton',...
                'String', 'Set Manual',...
                'Callback', @app.setManualCallback,...
                'Units', 'normalized',...
                'Position', [0.63, 0.525, 0.34, 0.45],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.peak1Txt = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'text',...
                'String', 'Peak 1:',...
                'Units', 'normalized',...
                'Position', [0, 0.025, 0.19, 0.3],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.peak1Edit = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'edit',...
                'Units', 'normalized',...
                'Position', [0.19, 0.1, 0.13, 0.25],...
                'Enable', 'off',...
                'Callback', @app.peak1Callback,...
                'FontSize', app.misc.fontSize);
            app.controls.peak2Txt = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'text',...
                'String', 'Peak 2:',...
                'Units', 'normalized',...
                'Position', [0.32, 0.025, 0.19, 0.3],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.peak2Edit = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'edit',...
                'Units', 'normalized',...
                'Position', [0.51, 0.1, 0.13, 0.25],...
                'Enable', 'off',...
                'Callback', @app.peak2Callback,...
                'FontSize', app.misc.fontSize);
            app.controls.minTxt = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'text',...
                'String', 'Min:',...
                'Units', 'normalized',...
                'Position', [0.64, 0.025, 0.19, 0.3],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.minEdit = uicontrol(app.controls.peakDetectMethodBtnGrp,...
                'Style', 'edit',...
                'Units', 'normalized',...
                'Position', [0.83, 0.1, 0.13, 0.25],...
                'Enable', 'off',...
                'Callback', @app.minCallback,...
                'FontSize', app.misc.fontSize);
            app.controls.whichBeadBtnGrp = uibuttongroup(app.controls.panel,...
                'Position', [0.01, 0.005, 0.48, 0.095],...
                'SelectionChangedFcn', @app.whichBeadCallback);
            app.controls.beadABtn = uicontrol(app.controls.whichBeadBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Bead A',...
                'Units', 'normalized',...
                'Position', [0.05, 0.05, 0.3, 0.9],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.beadBBtn = uicontrol(app.controls.whichBeadBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Bead B',...
                'Units', 'normalized',...
                'Position', [0.4, 0.05, 0.3, 0.9],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.bothBeadsBtn = uicontrol(app.controls.whichBeadBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Both',...
                'Value', 1,...
                'Units', 'normalized',...
                'Position', [0.75, 0.05, 0.2, 0.9],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            
            % Deselect panel (pop-up) of controls panel of tab 1:
            app.controls.deselect.panel = uipanel(app.controls.panel,...
                'Position', [0.02, 0.71, 0.46, 0.14],...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.deselect.panel.ForegroundColor = app.controls.deselect.panel.BackgroundColor;
            app.controls.deselect.deselectBtn = uicontrol(app.controls.panel,...
                'Style', 'togglebutton',...
                'String', 'Select Data for Removal',...
                'Callback', @app.deselectCallback,...
                'Units', 'normalized',...
                'Position', [0.02, 0.9, 0.46, 0.1],...
                'FontSize', app.misc.fontSize,...
                'Enable', 'off');
            app.controls.deselect.linePlot = plot(app.load.main.mainAxes, 0, 0, '.-',...
                'Color', app.misc.colors.deselect,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.controls.deselect.leftPlot = plot(app.load.main.mainAxes, 0, 0,...
                'LineWidth', 2,...
                'Color', 'none',...
                'HitTest', 'off');
            app.controls.deselect.rightPlot = plot(app.load.main.mainAxes, 0, 0,...
                'LineWidth', 2,...
                'Color', 'none',...
                'HitTest', 'off');
            app.controls.deselect.xpos = [];
            app.controls.deselect.patches = gobjects(0);
            app.controls.deselect.first = true;
            app.controls.deselect.undoBtn = uicontrol(app.controls.deselect.panel,...
                'Style', 'pushbutton',...
                'String', 'Undo Selection',...
                'Callback', @app.undoCallback,...
                'Units', 'normalized',...
                'Position', [0.025, 0.1, 0.35, 0.8],...
                'FontSize', app.misc.fontSize,...
                'Enable', 'off');
            app.controls.deselect.removeBtn = uicontrol(app.controls.deselect.panel,...
                'Style', 'pushbutton',...
                'String', 'Remove Data',...
                'Callback', @app.removeCallback,...
                'Units', 'normalized',...
                'Position', [0.4, 0.1, 0.35, 0.8],...
                'FontSize', app.misc.fontSize,...
                'Enable', 'off');
            app.controls.deselect.indices = [];
            app.controls.deselect.saveBtn = uicontrol(app.controls.deselect.panel,...
                'Style', 'pushbutton',...
                'String', 'Save',...
                'Callback', @app.saveCallback,...
                'Units', 'normalized',...
                'Position', [0.775, 0.1, 0.2, 0.8],...
                'FontSize', app.misc.fontSize,...
                'Enable', 'off');
            
            % Windows portion of controls panel of tab 1:
            app.controls.windows.panel = uipanel(app.misc.tab1,...
                'Position', [0.22, 0.03, 0.37, 0.45],...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.windows.windowsAxes = axes(app.controls.windows.panel,...
                'XLim', [1 2000],...
                'YLim', [3 2000],...
                'XGrid', 'on',...
                'YGrid', 'on',...
                'Box', 'on',...
                'NextPlot', 'add',...
                'YAxisLocation', 'right');
            xlabel(app.controls.windows.windowsAxes, 'Moving average filter window');
            ylabel(app.controls.windows.windowsAxes, '2nd order filter window');
            app.controls.windows.covwindow = 175; % In the paper, this is referred to as w_c.
            app.controls.windows.covsmooth = 73; % And this as w_s.
            app.controls.windows.defcovwindow = app.controls.windows.covwindow;
            app.controls.windows.defcovsmooth = app.controls.windows.covsmooth;
            app.controls.windows.covwindowTxt = uicontrol(app.controls.panel,...
                'Style', 'text',...
                'String', {'Moving average filter window:'; '(for calculating covariance)'},...
                'Units', 'normalized',...
                'Position', [0.5, 0.91, 0.38, 0.08],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.windows.covwindowEdit = uicontrol(app.controls.panel,...
                'Style', 'edit',...
                'String', app.controls.windows.defcovwindow,...
                'Units', 'normalized',...
                'Position', [0.88, 0.93, 0.1, 0.05],...
                'Enable', 'off',...
                'Callback', @app.covwindowCallback,...
                'FontSize', app.misc.fontSize);
            app.controls.windows.covsmoothTxt = uicontrol(app.controls.panel,...
                'Style', 'text',...
                'String', {'2nd order Savitzky-Golay filter window:'; '(for smoothing covariance)'},...
                'Units', 'normalized',...
                'Position', [0.5, 0.83, 0.38, 0.08],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.windows.covsmoothEdit = uicontrol(app.controls.panel,...
                'Style', 'edit',...
                'String', app.controls.windows.defcovsmooth,...
                'Units', 'normalized',...
                'Position', [0.88, 0.85, 0.1, 0.05],...
                'Enable', 'off',...
                'Callback', @app.covsmoothCallback,...
                'FontSize', app.misc.fontSize);
            app.controls.windows.interactBtn = uicontrol(app.controls.panel,...
                'Style', 'pushbutton',...
                'String', 'Choose Interactively',...
                'Callback', @app.windowsCallback,...
                'Units', 'normalized',...
                'Position', [0.53, 0.76, 0.2, 0.06],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.windows.pt = plot(app.controls.windows.windowsAxes, 0, 0, 'o',...
                'Color', app.misc.colors.windowsPt,...
                'MarkerSize', 10,...
                'Visible', 'off');
            app.controls.windows.pt.MarkerFaceColor = app.controls.windows.pt.Color;
            app.controls.windows.resetBtn = uicontrol(app.controls.panel,...
                'Style', 'pushbutton',...
                'String', 'Reset to Default',...
                'Callback', @app.defwindowsCallback,...
                'Units', 'normalized',...
                'Position', [0.76, 0.76, 0.2, 0.06],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            
            % Cutoffs portion of controls panel of tab 1:
            app.controls.cutoffs.defMinSep = app.minSep;
            app.controls.cutoffs.defMinDur = app.minDur;
            app.controls.cutoffs.minDurTxt = uicontrol(app.controls.panel,...
                'Style', 'text',...
                'String', {'Minimum event duration:'; '(for auto event removal)'},...
                'Units', 'normalized',...
                'Position', [0.5, 0.62, 0.25, 0.08],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.cutoffs.minDurEdit = uicontrol(app.controls.panel,...
                'Style', 'edit',...
                'String', app.minDur,...
                'Units', 'normalized',...
                'Position', [0.75, 0.64, 0.1, 0.05],...
                'Enable', 'off',...
                'Callback', @app.minDurCallback,...
                'FontSize', app.misc.fontSize);
            app.controls.cutoffs.minSepTxt = uicontrol(app.controls.panel,...
                'Style', 'text',...
                'String', {'Minimum event separation:'; '(for auto event removal)'},...
                'Units', 'normalized',...
                'Position', [0.5, 0.54, 0.25, 0.08],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            app.controls.cutoffs.minSepEdit = uicontrol(app.controls.panel,...
                'Style', 'edit',...
                'String', app.minSep,...
                'Units', 'normalized',...
                'Position', [0.75, 0.56, 0.1, 0.05],...
                'Enable', 'off',...
                'Callback', @app.minSepCallback,...
                'FontSize', app.misc.fontSize);
            app.controls.cutoffs.resetBtn = uicontrol(app.controls.panel,...
                'Style', 'pushbutton',...
                'String', '<html><p align="center">Reset to<br/>Default</p></html>',...
                'Callback', @app.defcutoffsCallback,...
                'Units', 'normalized',...
                'Position', [0.87, 0.565, 0.11, 0.12],...
                'Enable', 'off',...
                'FontSize', app.misc.fontSize);
            
            % allEvents becomes outdated when the covariance/peaks/min/bead(s) change.
            % Updating allEvents requires running the changepoint algorithm, which is
            % expensive, so the user is always asked before it happens. Therefore, the
            % values of the covariance/peaks/min/whichBeads that were used to determine
            % allEvents need to be stored so current values can be compared.
            app.misc.active.allEvents.cov = [];
            app.misc.active.allEvents.peak1 = [];
            app.misc.active.allEvents.peak2 = [];
            app.misc.active.allEvents.min = [];
            app.misc.active.allEvents.whichBeads = app.controls.whichBeadBtnGrp.SelectedObject;
            app.misc.active.allEvents.useChangepoint = app.controls.useChangepointBtn.Value;
            
            % Save to excel panel of tab 1:
            app.excel.panel = uipanel(app.controls.panel,...
                'Title', 'Save to Excel',...
                'Position', [0.5, 0.005, 0.5, 0.5],...
                'FontSize', app.misc.fontSize);
            app.excel.listbox = uicontrol(app.excel.panel,...
                'Style', 'listbox',...
                'String', strings,...
                'Value', [],...
                'Units', 'normalized',...
                'Position', [0.05, 0.2, 0.9, 0.77],...
                'Max', 2,... This allows user to select multiple files using ctrl+click.
                'FontSize', app.misc.fontSize);
            app.excel.files = strings;
            app.excel.createBtn = uicontrol(app.excel.panel,...
                'Style', 'pushbutton',...
                'String', 'New',...
                'Callback', @app.createExcelCallback,...
                'Units', 'normalized',...
                'Position', [0.05, 0.02, 0.25, 0.15],...
                'FontSize', app.misc.fontSize);
            app.excel.findBtn = uicontrol(app.excel.panel,...
                'Style', 'pushbutton',...
                'String', 'Find',...
                'Callback', @app.findExcelCallback,...
                'Units', 'normalized',...
                'Position', [0.33, 0.02, 0.25, 0.15],...
                'FontSize', app.misc.fontSize);
            app.excel.addBtn = uicontrol(app.excel.panel,...
                'Style', 'pushbutton',...
                'String', 'Save',...
                'Callback', @app.addToExcelCallback,...
                'Enable', 'off',...
                'Units', 'normalized',...
                'Position', [0.61, 0.02, 0.34, 0.15],...
                'FontSize', app.misc.fontSize);
            [app.excel.cell.path,...
                app.excel.cell.fs,...
                app.excel.cell.N,...
                app.excel.cell.flipped,...
                app.excel.cell.deselectIndices,...
                app.excel.cell.minDur,...
                app.excel.cell.minSep,...
                app.excel.cell.covwindow,...
                app.excel.cell.covsmooth,...
                app.excel.cell.autopeak1,...
                app.excel.cell.autopeak2,...
                app.excel.cell.automin,...
                app.excel.cell.manpeak1,...
                app.excel.cell.manpeak2,...
                app.excel.cell.manmin,...
                app.excel.cell.usemin,...
                app.excel.cell.useman,...
                app.excel.cell.useChngpt,...
                app.excel.cell.whichBeads,...
                app.excel.cell.colNames,...
                app.excel.cell.data] = app.dealCells;
            
            % Load panel of tab 2:
            app.load2.panel = uipanel(app.misc.tab2,...
                'Title', 'Load',...
                'Position', [0.005, 0.52, 0.99, 0.48],...
                'FontSize', app.misc.fontSize);
            app.load2.inputBtn = uicontrol(app.load2.panel,...
                'Style', 'pushbutton',...
                'String', 'Choose Files',...
                'Units', 'normalized',...
                'Position', [0.005, 0.902, 0.085, 0.083],...
                'FontSize', app.misc.fontSize,...
                'Callback', @app.chooseXlsxCallback);
            app.load2.table = uitable(app.load2.panel,...
                'Units', 'normalized',...
                'RowName', [],...
                'ColumnName', {'Sheet', 'File', '# Events', 'Include'},...
                'ColumnFormat', [repmat({'char'}, 1, 3), 'logical'],...
                'ColumnEditable', [false(1, 3), true],...
                'Position', [0.01, 0.05, 0.98, 0.81],...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load2.tableData = [];
            app.load2.maxHeight = app.load2.table.Position(4);
            app.load2.tWidths = [];
            app.data2.fs = [];
            app.data2.events = [];
            app.data2.beadA = [];
            app.data2.beadB = [];
            app.load2.whichBeadBtnGrp = uibuttongroup(app.load2.panel,...
                'Position', [0.095, 0.902, 0.2, 0.083],...
                'SelectionChangedFcn', @app.whichBeadCallback2,...
                'Visible', 'off');
            app.load2.beadABtn = uicontrol(app.load2.whichBeadBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Bead A',...
                'Units', 'normalized',...
                'Position', [0.05, 0.05, 0.3, 0.9],...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load2.beadBBtn = uicontrol(app.load2.whichBeadBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Bead B',...
                'Units', 'normalized',...
                'Position', [0.4, 0.05, 0.3, 0.9],...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load2.bothBeadsBtn = uicontrol(app.load2.whichBeadBtnGrp,...
                'Style', 'radiobutton',...
                'String', 'Both',...
                'Value', 1,...
                'Units', 'normalized',...
                'Position', [0.75, 0.05, 0.2, 0.9],...
                'Visible', 'off',...
                'FontSize', app.misc.fontSize);
            app.load2.currentTask = uicontrol(app.load2.panel,...
                'Style', 'text',...
                'String', '',...
                'ForegroundColor', app.misc.colors.currentTask,...
                'Units', 'normalized',...
                'Position', [0.99, 0.93, 0.01, 0.055],...
                'FontSize', 1.3*app.misc.fontSize,...
                'HorizontalAlignment', 'right');
            app.load2.taskMaxWidth = 0.61;
            app.load2.skipBtn = uicontrol(app.load2.panel,...
                'Style', 'pushbutton',...
                'String', 'Skip',...
                'Units', 'normalized',...
                'Position', [0.3, 0.902, 0.085, 0.083],...
                'FontSize', app.misc.fontSize,...
                'Callback', @app.skipBtnCallback,...
                'Visible', 'off');
            
            % Analyze panel of tab 2:
            app.analyze2.panel = uipanel(app.misc.tab2,...
                'Title', 'Analyze',...
                'Position', [0.005, 0.005, 0.99, 0.513],...
                'FontSize', app.misc.fontSize);
            
            % Miscellaneous axes of analyze panel of tab 2:
            app.analyze2.misc.miscAxes = axes(app.analyze2.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.01, 0.092, 0.45, 0.89],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'combinedMiscAxes');
            app.analyze2.misc.miscPlot = plot(app.analyze2.misc.miscAxes, 0, 0,...
                'Color', app.misc.colors.misc,...
                'MarkerSize', 3,...
                'MarkerFaceColor', app.misc.colors.misc,...
                'MarkerEdgeColor', app.misc.colors.misc,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.misc.miscFitPlot = plot(app.analyze2.misc.miscAxes, 0, 0,...
                'Color', app.misc.colors.miscFit,...
                'LineWidth', 1.5,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.misc.miscLabel = text(app.analyze2.misc.miscAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.misc.data = 'Event Duration'; % Show event durations by default.
            app.data2.durs = [];
            app.data2.beadAStep1 = [];
            app.data2.beadBStep1 = [];
            app.data2.beadAStep2 = [];
            app.data2.beadBStep2 = [];
            app.data2.beadATotalStep = [];
            app.data2.beadBTotalStep = [];
            app.data2.forceA = [];
            app.data2.forceB = [];
            app.analyze2.misc.kDur = [];
            app.analyze2.misc.meanStep1 = [];
            app.analyze2.misc.stdStep1 = [];
            app.analyze2.misc.meanStep2 = [];
            app.analyze2.misc.stdStep2 = [];
            app.analyze2.misc.meanTotalStep = [];
            app.analyze2.misc.stdTotalStep = [];
            
            % Ensemble average axes of analyze panel of tab 2:
            app.analyze2.ensemble.ensembleAxes = axes(app.analyze2.panel,...
                'Box', 'on',...
                'Visible', 'off',...
                'Position', [0.51, 0.092, 0.45, 0.89],...
                'FontSize', app.misc.fontSize,...
                'NextPlot', 'add',...
                'YAxisLocation', 'right',...
                'UserData', 'combinedEnsembleAxes');
            xlabel(app.analyze2.ensemble.ensembleAxes, 'time (s)');
            app.analyze2.ensemble.globalFit = true;
            app.analyze2.ensemble.startPlot = plot(app.analyze2.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensStart,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.ensemble.endPlot = plot(app.analyze2.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensEnd,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.ensemble.startFitPlot = plot(app.analyze2.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensStartFit,...
                'LineWidth', 2,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.ensemble.endFitPlot = plot(app.analyze2.ensemble.ensembleAxes, 0, 0,...
                'Color', app.misc.colors.ensEndFit,...
                'LineWidth', 2,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.ensemble.startLabel = text(app.analyze2.ensemble.ensembleAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.ensemble.endLabel = text(app.analyze2.ensemble.ensembleAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.ensemble.stepLabel = text(app.analyze2.ensemble.ensembleAxes, 0, 0, '',...
                'FontSize', app.misc.fontSize,...
                'Visible', 'off',...
                'HitTest', 'off');
            app.analyze2.ensemble.start = [];
            app.analyze2.ensemble.end = [];
            app.analyze2.ensemble.startFit = [];
            app.analyze2.ensemble.endFit = [];
            app.analyze2.ensemble.p = [];
            
            if isprop(app.load.main.mainAxes, 'Toolbar')
                % R2018b or later
                app.misc.axToolbars = true;
                app.setupAxToolbar(app.load.main.mainAxes);
                app.setupAxToolbar(app.load.cov.covAxes);
                app.setupAxToolbar(app.analyze.covHist.covHistAxes);
                app.setupAxToolbar(app.analyze.misc.miscAxes);
                app.setupAxToolbar(app.analyze.ensemble.ensembleAxes);
                app.setupAxToolbar(app.analyze.individual.likelihoodAxes);
                app.setupAxToolbar(app.analyze.individual.eventAxes);
                app.setupAxToolbar(app.analyze2.misc.miscAxes);
                app.setupAxToolbar(app.analyze2.ensemble.ensembleAxes);
            else
                % pre-R2018b
                app.misc.axToolbars = false;
                app.misc.fig.MenuBar = 'figure';
                menus = findall(app.misc.fig, 'type', 'uimenu');
                delete(findobj(menus, '-not', 'tag', 'figMenuZoomIn', '-and',...
                    '-not', 'tag', 'figMenuZoomOut', '-and',...
                    '-not', 'tag', 'figMenuPan', '-and',...
                    '-not', 'tag', 'figMenuRotate3D', '-and',...
                    '-not', 'tag', 'figMenuResetView', '-and',...
                    '-not', 'tag', 'figMenuFileSaveAs', '-and',...
                    '-not', 'tag', 'figMenuFile', '-and',...
                    '-not', 'tag', 'figMenuTools'));
                app.misc.submenus = findall(app.misc.fig, 'type', 'uimenu', '-and',...
                    '-not', 'tag', 'figMenuTools', '-and',...
                    '-not', 'tag', 'figMenuFile');
                app.misc.menus = findall(app.misc.fig, 'type', 'uimenu', '-and',...
                    'tag', 'figMenuTools', '-or', 'tag', 'figMenuFile');
                set(app.misc.submenus, 'MenuSelectedFcn', @app.menuSelectedFcn);
            end
            
            % A double array where each index corresponds to an 'enableable' uicontrol, or
            % a uicontrol which is disabled while the program is busy. 1 in this array
            % indicates the corresponding uicontrol should become enabled after the
            % program finishes. -1 indicates the uicontrol should become disabled. 0 means
            % the uicontrol will return to whatever state it was in before the program
            % became busy.
            app.misc.toBeEnabled = zeros(length(app.getEnableable), 1);
            app.misc.disabled = false;
            
            % This determines the properties for which SetObservable is true...
            c = ?SPASM;
            observable = cell(1, length(c.PropertyList));
            for i = 1:length(c.PropertyList)
                if c.PropertyList(i).SetObservable
                    observable{i} = c.PropertyList(i).Name;
                end
            end
            observable = observable(~cellfun(@isempty, observable));
            % ...and attaches a listener which will call postSetCallback() after any of
            % the properties change value.
            app.misc.postListener = addlistener(app, observable, 'PostSet', @app.postSetCallback);
            app.misc.postListener.Recursive = true;
            app.misc.listenerError = [];
            
            app.misc.signal = app.checkInstall('Signal Processing Toolbox');
            app.misc.optim = app.checkInstall('Optimization Toolbox');
            
            % Final touches.
            movegui(app.misc.fig, 'center')
            app.misc.fig.Visible = 'on';
            app.misc.fig.SizeChangedFcn = @app.sizeChange;
            
            if ~app.misc.signal
                warndlg('The Signal Processing Toolbox is needed to smooth the covariance.')
            end
            
            if ~app.misc.optim
                warndlg('The Optimization Toolbox is needed to determine fits for the ensemble averages.')
            end
        end
        
        function closeRequest(app, ~, ~)
            % Makes sure to delete the app if user closes the figure.
            
            closereq;
            delete(app);
        end
        
        function sizeChange(app, ~, ~)
            % The default sizeChanged callback for the figure. Most elements resize
            % automatically, but this function updates the font sizes as needed.
            
            try
                oldFS = app.misc.fontSize;
                app.misc.fontSize = min(app.misc.fig.Position(4) / 101.25, app.misc.fig.Position(3) / 180);
                set(findall(app.misc.fig, 'FontSize', oldFS), 'FontSize', app.misc.fontSize);
                app.load.currentTask.FontSize = 1.3*app.misc.fontSize;
                app.load2.currentTask.FontSize = 1.3*app.misc.fontSize;
                if strcmp(app.load.KATxt.Visible, 'on')
                    app.applyHeader; % Update labels whose positions depend on their sizes.
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.SizeChangedFcn = [];
                app.unexpectedError(ME, sprintf(['An unexpected error occurred while resizing the figure:\n\n%s\n\n'...
                    'Exiting the app...'], ME.message));
                close(app.misc.fig);
                delete(app);
            end
        end
        
        function defaultMotion(app, ~, ~)
            % The default motion callback for the figure. It displays a label when the
            % user mouses over event patches.
            
            try 
                if app.showEvents
                    pi = app.load.main.prevIn;
                    if ~isempty(pi) && ~isequal(app.getCurrentEventAsIndex, pi)
                        % If mouse was previously in an inactive event which is not the
                        % current event, the patch should be hidden.
                        app.load.main.patches(pi).Visible = 'off';
                        app.load.main.prevIn = [];
                    end
                    t = app.load.main.label; % Text label.
                    [in, p1] = app.inAxes(app.load.main.mainAxes); % [whether in axes, mouse position in axes]
                    if in && isequal(app.misc.tab1, app.misc.tabgroup.SelectedTab)
                        fp = get(app.misc.fig, 'CurrentPoint'); % Mouse position in pixels.
                        pp = getpixelposition(app.load.panel); % Load panel position in tab 1.
                        up = getpixelposition(app.misc.tab1); % Tab 1 position in figure.
                        fx = fp(1); % Mouse X in figure.
                        fy = fp(2); % Mouse Y in figure.
                        ps = app.load.main.patches; % Event patches.
                        j = 0; k = 0; % To determine the index of the event patch that mouse is currently in.
                        for i = 1:length(ps)
                            % For each patch...
                            try
                                p = ps(i);
                                if strcmp(p.Visible, 'off') && ~app.showRemoved
                                    % ...if patch is invisible (i.e. removed) and removed
                                    % events should be invisible, skip.
                                    continue
                                end
                                % ...otherwise, increment j.
                                j = j + 1;
                                if inpolygon(p1(1), p1(2), p.XData, p.YData)
                                    % If mouse is in this patch, store value of j in k.
                                    k = j;
                                    break
                                end
                            catch
                            end
                        end
                        if k ~= 0
                            % If mouse is in some patch...
                            switch app.showRemoved
                                case 1
                                    % ...removed events should be shown, so show this patch.
                                    % Also update the label text.
                                    ps(k).Visible = 'on';
                                    t.String = [num2str(k) '/' num2str(length(app.activeEvents))];
                                    if ~app.activeEvents(k)
                                        % If this patch is removed, remember its index so that
                                        % it will not stay visible should mouse leave.
                                        app.load.main.prevIn = k;
                                    end
                                case 0
                                    % ...removed events are not being shown, so this patch
                                    % must already be showing. Just update the label text.
                                    t.String = [num2str(k) '/' num2str(sum(app.activeEvents))];
                            end
                            t.Position(1:2) = [fx-pp(1)-up(1), fy-pp(2)-up(2)]; % Update the label position.
                            t.FitBoxToText = 'on';
                            t.Visible = 'on';
                        else
                            % Mouse is not in any patches.
                            t.Visible = 'off';
                        end
                    else
                        % Mouse is not even in axes.
                        t.Visible = 'off';
                    end
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowButtonMotionFcn = [];
                app.unexpectedError(ME, sprintf(['An unexpected error occurred while handling cursor movement:\n\n%s\n\n'...
                    'Exiting the app...'], ME.message));
                close(app.misc.fig);
                delete(app);
            end
        end
        
        function defaultClick(app, ~, ~)
            % The default button down callback for the figure. It determines if the user
            % clicked in any axes and, if so, gets the UIContextMenu ready for the axes
            % which were clicked.
            
            try
                if app.misc.disabled
                    enableState = 'off';
                else
                    enableState = 'on';
                end

                switch app.misc.tabgroup.SelectedTab
                    % Only check axes which are on the tab currently being shown.
                    case app.misc.tab1
                        axes = [app.load.main.mainAxes,...
                            app.load.cov.covAxes,...
                            app.analyze.covHist.covHistAxes,...
                            app.analyze.misc.miscAxes,...
                            app.analyze.ensemble.ensembleAxes,...
                            app.analyze.individual.likelihoodAxes,...
                            app.analyze.individual.eventAxes];
                    case app.misc.tab2
                        axes = [app.analyze2.misc.miscAxes,...
                            app.analyze2.ensemble.ensembleAxes];
                end

                % Determine which sets of axes the mouse is inside.
                in = app.inAxes(axes);
                ax = axes(in);
                a = ax(arrayfun(@(a) strcmp(a.Visible,'on'), ax)); % Only consider axes which are visible.

                % Set the UIContextMenu.
                if ~isempty(a)

                    % Just in case there are multiple axes in a, choose the first. (But there
                    % shouldn't be.)
                    a = a(1);

                    % Delete all previous submenus of the UIContextMenu.
                    children = app.misc.cm.Children;
                    for c = 1:numel(children)
                        delete(children(c)); 
                    end

                    % Add a submenu for saving the figure as a .fig file, passing the clicked
                    % axes a into the submenu's callback.
                    uimenu(app.misc.cm,...
                        'Text', 'Save as .fig',...
                        'MenuSelectedFcn', {@app.saveFigCallback, a},...
                        'Enable', enableState);

                    % Add a submenu for saving the figure as a .eps file, passing the clicked
                    % axes a into the submenu's callback.
                    uimenu(app.misc.cm,...
                        'Text', 'Save as .eps',...
                        'MenuSelectedFcn', {@app.saveEpsCallback, a},...
                        'Enable', enableState);

                    % If user clicked a set of axes showing ensemble averages, also add a
                    % submenu for exporting the data to Excel.
                    if strcmp(a.UserData, 'ensembleAxes') || strcmp(a.UserData, 'combinedEnsembleAxes')
                        uimenu(app.misc.cm,...
                            'Text', 'Export Data',...
                            'MenuSelectedFcn', {@app.exportEnsembleCallback, a},...
                            'Enable', enableState);
                        if strcmp(a.UserData, 'ensembleAxes')
                            globalFit = app.analyze.ensemble.globalFit;
                        else
                            globalFit = app.analyze2.ensemble.globalFit;
                        end
                        checked = 'off';
                        if globalFit
                            checked = 'on';
                        end
                        uimenu(app.misc.cm,...
                            'Text', 'Fit Averages Globally',...
                            'Checked', checked,...
                            'MenuSelectedFcn', {@app.fitGloballyCallback, a},...
                            'Enable', enableState);
                    end

                    % If user clicked one of the misc axes, also add a submenu letting them
                    % choose which data to show in the plot.
                    if strcmp(a.UserData, 'miscAxes') || strcmp(a.UserData, 'combinedMiscAxes')
                        uimenu(app.misc.cm,...
                            'Text', 'Export Data',...
                            'MenuSelectedFcn', {@app.exportMiscDistCallback, a},...
                            'Enable', enableState);
                        if strcmp(a.UserData, 'miscAxes')
                            s = app.analyze.misc.data;
                            disableForceDur = isequal(app.misc.active.allEvents.whichBeads, app.controls.bothBeadsBtn);
                        else
                            s = app.analyze2.misc.data;
                            disableForceDur = app.load2.bothBeadsBtn.Value;
                        end
                        m = uimenu(app.misc.cm,...
                            'Text', 'Choose Data',...
                            'Enable', enableState);
                        addMenu('Event Duration');
                        addMenu('Step 1 Size');
                        addMenu('Step 2 Size');
                        addMenu('Total Step Size');
                        addMenu('Event Duration / Force on Selected Bead', disableForceDur);
                    end

                    % Attach the UIContextMenu to a.
                    a.UIContextMenu = app.misc.cm;
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowButtonDownFcn = [];
                app.unexpectedError(ME, sprintf(['An unexpected error occurred while handling a mouse click:\n\n%s\n\n'...
                    'Exiting the app...'], ME.message));
                close(app.misc.fig);
                delete(app);
            end
            
            function addMenu(str, disable)
                checked = 'on';
                if ~strcmp(s, str)
                    checked = 'off';
                end
                if nargin == 2 && disable
                    enable = false;
                else
                    enable = enableState;
                end
                uimenu(m,...
                    'Text', str,...
                    'MenuSelectedFcn', {@app.chooseDataCallback, a},...
                    'Checked', checked,...
                    'Enable', enable);
            end
        end
        
        function saveFigCallback(app, ~, ~, ax)
            % Called when user clicks 'Save Figure' button on the UIContextMenu of ax.
            % Opens a dialog box allowing them to save a figure containing ax as a .fig
            % file.
            
            try
                if app.misc.disabled
                    return
                end

                store = app.disable;
                drawnow;

                [name, path] = uiputfile([ax.UserData '.fig']);
                if name == 0
                    app.restore(store);
                    return
                end
                [~, name] = fileparts([path name]);
                app.updateCurrentTask(sprintf('Saving %s.fig.', name));

                % Create a temporary figure whose dimensions are slightly bigger than ax, and
                % copy ax onto that figure.
                f = figure('Visible', 'off');
                ax.Units = 'pixels';
                f.Position(3) = ax.Position(3) + ax.TightInset(1) + ax.TightInset(3) + 10;
                f.Position(4) = ax.Position(4) + ax.TightInset(2) + ax.TightInset(4) + 10;
                h = copyobj(ax, f);
                h.Position([1 2]) = [5+ax.TightInset(1) 5+ax.TightInset(2)];
                ax.Units = 'normalized';
                children = h.Children;
                for c = 1:numel(children)
                    % Turn HitTest of children back on, so user can select all plot elements.
                    if isprop(children(c), 'HitTest')
                        children(c).HitTest = 'on';
                    end
                    % Delete any invisible children.
                    if isprop(children(c), 'Visible') && strcmp(children(c).Visible, 'off')
                        delete(children(c));
                    end
                end

                % The temporary figure is invisible so that user cannot see it. Make sure that
                % the saved figure will be visible when the .fig file is opened.
                f.CreateFcn = @(src, ~) set(src, 'Visible', 'on');
                
                % Save the .fig file.
                savefig(f, [path name '.fig']);
                
                % Delete the temporary figure.
                delete(f);
                
                app.restore(store);
            catch ME
                try
                    delete(f);
                catch
                end
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME, sprintf('An unexpected error occurred while saving:\n\n%s', ME.message));
                app.restore(store);
            end
        end
        
        function saveEpsCallback(app, ~, ~, ax)
            % Called when user clicks 'Save Figure' button on the UIContextMenu of ax.
            % Opens a dialog box allowing them to save a figure containing ax as a .eps
            % file.
            
            try
                if app.misc.disabled
                    return
                end

                store = app.disable;
                drawnow;

                [name, path] = uiputfile([ax.UserData '.eps']);
                if name == 0
                    app.restore(store);
                    return
                end
                [~, name] = fileparts([path name]);
                app.updateCurrentTask(sprintf('Saving %s.eps.', name));

                % Create a temporary figure whose dimensions are slightly bigger than ax, and
                % copy ax onto that figure.
                f = figure('Visible', 'off');
                ax.Units = 'pixels';
                f.Position(3) = ax.Position(3) + ax.TightInset(1) + ax.TightInset(3) + 10;
                f.Position(4) = ax.Position(4) + ax.TightInset(2) + ax.TightInset(4) + 10;
                h = copyobj(ax, f);
                h.Position([1 2]) = [5+ax.TightInset(1) 5+ax.TightInset(2)];
                ax.Units = 'normalized';
                children = h.Children;
                for c = 1:numel(children)
                    % Turn HitTest of children back on, so user can select all plot elements.
                    if isprop(children(c), 'HitTest')
                        children(c).HitTest = 'on';
                    end
                    % Delete any invisible children.
                    if isprop(children(c), 'Visible') && strcmp(children(c).Visible, 'off')
                        delete(children(c));
                    end
                end

                % The temporary figure is invisible so that user cannot see it. Make sure that
                % the saved figure will be visible when the .fig file is opened.
                f.CreateFcn = @(src, ~) set(src, 'Visible', 'on');

                % Save the .fig file.
                saveas(f, [path name], 'epsc');

                % Delete the temporary figure.
                delete(f);

                app.restore(store);
            catch ME
                try
                    delete(f);
                catch
                end
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME, sprintf('An unexpected error occurred while saving:\n\n%s', ME.message));
                app.restore(store);
            end
        end
        
        function exportEnsembleCallback(app, ~, ~, ax)
            % Called when user clicks 'Export Data' button on the UIContextMenu of ax.
            % Opens a dialog box allowing them to save the corresponding ensemble average
            % data to a .xlsx file.
            
            try
                if app.misc.disabled
                    return
                end

                store = app.disable;
                drawnow;

                [name, path] = uiputfile([ax.UserData 'Data.xlsx']);
                if name == 0
                    app.restore(store);
                    return
                end

                % Obtain the filename without the extension, in case user changed the
                % extension in the dialog box.
                [~, name] = fileparts([path name]);
                app.updateCurrentTask(sprintf('Saving %s.xlsx.', name));

                % If user picks a file which already exists, first delete the file (because
                % writetable will not overwrite, at least not in R2019a).
                if exist([path name '.xlsx'], 'file') == 2
                  delete([path name '.xlsx']);
                end
                
                % Determine the data which corresponds to the selected axes.
                switch ax.UserData
                    case 'combinedEnsembleAxes'
                        fs = app.data2.fs;
                        pts_bef = round(app.analyze.ensemble.pts_bef * fs);
                        pts_after = round(app.analyze.ensemble.pts_after * fs);
                        skip = round(app.analyze.ensemble.num_pts_to_skip);
                        
                        startAvg = app.analyze2.ensemble.startPlot.YData;
                        startTime = ((1:numel(startAvg)) - pts_bef - 1) / fs;
                        
                        endAvg = app.analyze2.ensemble.endPlot.YData;
                        endTime = ((-numel(endAvg):-1) + pts_after) / fs;
                        
                        startFit = app.analyze2.ensemble.startFitPlot.YData;
                        startFitTime = (skip + (1:numel(startFit))) / fs;
                        if isequal(startFit, [0, 0]) || strcmp(app.analyze2.ensemble.startFitPlot.Visible, 'off')
                            startFit = NaN;
                        end
                        
                        endFit = app.analyze2.ensemble.endFitPlot.YData;
                        endFitTime = ((-numel(endFit):-1) - skip) / fs;
                        if isequal(endFit, [0, 0]) || strcmp(app.analyze2.ensemble.endFitPlot.Visible, 'off')
                            endFit = NaN;
                        end
                        
                        if isempty(app.analyze2.ensemble.p)
                            kf = NaN;
                            kr = NaN;
                            step1 = NaN;
                            step2 = NaN;
                            totalStep = NaN;
                        else
                            kf = app.analyze2.ensemble.p(1);
                            kr = app.analyze2.ensemble.p(2);
                            step1 = app.analyze2.ensemble.p(3);
                            totalStep = app.analyze2.ensemble.p(4);
                            step2 = totalStep - step1;
                        end
                    otherwise
                        fs = app.data.fs;
                        pts_bef = round(app.analyze.ensemble.pts_bef * fs);
                        pts_after = round(app.analyze.ensemble.pts_after * fs);
                        skip = round(app.analyze.ensemble.num_pts_to_skip);
                        
                        startAvg = app.analyze.ensemble.startPlot.YData;
                        startTime = ((1:numel(startAvg)) - pts_bef - 1) / fs;
                        
                        endAvg = app.analyze.ensemble.endPlot.YData;
                        endTime = ((-numel(endAvg):-1) + pts_after) / fs;
                        
                        startFit = app.analyze.ensemble.startFitPlot.YData;
                        startFitTime = (skip + (1:numel(startFit))) / fs;
                        if isequal(startFit, [0, 0]) || strcmp(app.analyze.ensemble.startFitPlot.Visible, 'off')
                            startFit = NaN;
                            startFitTime = NaN;
                        end
                        
                        endFit = app.analyze.ensemble.endFitPlot.YData;
                        endFitTime = ((-numel(endFit):-1) - skip) / fs;
                        if isequal(endFit, [0, 0]) || strcmp(app.analyze.ensemble.endFitPlot.Visible, 'off')
                            endFit = NaN;
                            endFitTime = NaN;
                        end
                        
                        if isempty(app.analyze.ensemble.p)
                            kf = NaN;
                            kr = NaN;
                            step1 = NaN;
                            step2 = NaN;
                            totalStep = NaN;
                        else
                            kf = app.analyze.ensemble.p(1);
                            kr = app.analyze.ensemble.p(2);
                            step1 = app.analyze.ensemble.p(3);
                            totalStep = app.analyze.ensemble.p(4);
                            step2 = totalStep - step1;
                        end
                        if app.flipped
                            step1 = -step1;
                            step2 = -step2;
                            totalStep = -totalStep;
                        end
                end
                colNames = {'start average time (s)'...
                    'start average position (nm)'...
                    'start fit time (s)'...
                    'start fit position (nm)'...
                    ''...
                    'end average time (s)'...
                    'end average position (nm)'...
                    'end fit time (s)'...
                    'end fit position (s)'};
                rowNames = {'time-forward k (1/s)';
                    'time-reverse k (1/s)';
                    'step 1 (nm)';
                    'step 2 (nm)';
                    'total step (nm)'};
                
                % Write the data to the specified .xlsx file.
                writetable(table(colNames),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(1) '1']);
                writetable(table([startTime(:), startAvg(:)]),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(1) '2']);
                writetable(table([startFitTime(:), startFit(:)]),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(3) '2']);
                writetable(table([endTime(:), endAvg(:)]),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(6) '2']);
                writetable(table([endFitTime(:), endFit(:)]),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(8) '2']);
                writetable(table(rowNames),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(11) '1']);
                writetable(table([kf; kr; step1; step2; totalStep]),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(12) '1']);

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME, sprintf('An unexpected error occurred while saving:\n\n%s', ME.message));
                app.restore(store);
            end
        end
        
        function fitGloballyCallback(app, src, ~, ax)
            % Called when user toggles 'Fit Averages Globally' button on the UIContextMenu
            % of ax. Recalculates the ensemble average fits accordingly.
            
            switch ax.UserData
                case 'ensembleAxes'
                    backup = app.storeBackupTab1;
                case 'combinedEnsembleAxes'
                    backup = app.storeBackupTab2;
            end
            
            try
                if app.misc.disabled
                    return
                end

                store = app.disable;
                drawnow;
                
                % Determine the data which corresponds to the selected axes.
                switch ax.UserData
                    case 'combinedEnsembleAxes'
                        app.analyze2.ensemble.globalFit = ~app.analyze2.ensemble.globalFit;
                        if app.analyze2.ensemble.globalFit
                            src.Checked = 'on';
                        else
                            src.Checked = 'off';
                        end
                        if app.misc.optim
                            [app.analyze2.ensemble.startFit, app.analyze2.ensemble.endFit,...
                                app.analyze2.ensemble.p] = app.calcEnsembleFits2;
                            app.plotEnsembleFits2;
                        end
                    otherwise
                        app.analyze.ensemble.globalFit = ~app.analyze.ensemble.globalFit;
                        if app.analyze.ensemble.globalFit
                            src.Checked = 'on';
                        else
                            src.Checked = 'off';
                        end
                        if app.misc.optim
                            [app.analyze.ensemble.startFit, app.analyze.ensemble.endFit,...
                                app.analyze.ensemble.p] = app.calcEnsembleFits;
                            app.plotEnsembleFits;
                        end
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                switch ax.UserData
                    case 'ensembleAxes'
                        app.restoreBackupTab1(backup);
                    case 'combinedEnsembleAxes'
                        app.restoreBackupTab2(backup);
                end
            end
        end
        
        function exportMiscDistCallback(app, ~, ~, ax)
            % Called when user clicks 'Export Data' button on the UIContextMenu of ax.
            % Opens a dialog box allowing them to save the data and, if present, fit of
            % corresponding misc axes to a .xlsx file.
            
            try
                if app.misc.disabled
                    return
                end

                store = app.disable;
                drawnow;

                % Determine which plot is currently showing, and create suggested filename.
                switch ax.UserData
                    case 'combinedMiscAxes'
                        showing = app.analyze2.misc.data;
                        name = ['combined' strrep(strrep(showing,' ',''),'/','') 'Data.xlsx'];
                    otherwise
                        showing = app.analyze.misc.data;
                        name = [strrep(strrep(showing,' ',''),'/','') 'Data.xlsx'];
                end

                [name, path] = uiputfile(name);
                if name == 0
                    app.restore(store);
                    return
                end

                % Obtain the filename without the extension, in case user changed the
                % extension in the dialog box.
                [~, name] = fileparts([path name]);
                app.updateCurrentTask(sprintf('Saving %s.xlsx.', name));

                % If user picks a file which already exists, first delete the file (because
                % writetable will not overwrite, at least not in R2019a).
                if exist([path name '.xlsx'], 'file') == 2
                  delete([path name '.xlsx']);
                end

                % Determine the data which corresponds to the selected axes.
                switch ax.UserData
                    case 'combinedMiscAxes'
                        miscData = [app.analyze2.misc.miscPlot.XData(:),...
                            app.analyze2.misc.miscPlot.YData(:)];
                        miscFit = [app.analyze2.misc.miscFitPlot.XData(:),...
                            app.analyze2.misc.miscFitPlot.YData(:)];
                        kDur = abs(app.analyze2.misc.kDur);
                        meanStep1 = app.analyze2.misc.meanStep1;
                        stdStep1 = app.analyze2.misc.stdStep1;
                        meanStep2 = app.analyze2.misc.meanStep2;
                        stdStep2 = app.analyze2.misc.stdStep2;
                        meanTotalStep = app.analyze2.misc.meanTotalStep;
                        stdTotalStep = app.analyze2.misc.stdTotalStep;
                    otherwise
                        miscData = [app.analyze.misc.miscPlot.XData(:),...
                            app.analyze.misc.miscPlot.YData(:)];
                        miscFit = [app.analyze.misc.miscFitPlot.XData(:),...
                            app.analyze.misc.miscFitPlot.YData(:)];
                        kDur = abs(app.analyze.misc.kDur);
                        meanStep1 = app.analyze.misc.meanStep1;
                        stdStep1 = app.analyze.misc.stdStep1;
                        meanStep2 = app.analyze.misc.meanStep2;
                        stdStep2 = app.analyze.misc.stdStep2;
                        meanTotalStep = app.analyze.misc.meanTotalStep;
                        stdTotalStep = app.analyze.misc.stdTotalStep;
                end
                
                % Determine labels and parameters to output based on which plot is currently
                % showing.
                switch showing
                    case 'Step 1 Size'
                        T = table([meanStep1; stdStep1]);
                        colNames = {'step 1 sizes (nm)','probability'...
                            'fit sizes (nm)','fit probability'};
                        rowNames = {'mean (nm)';'standard deviation (nm)'};
                    case 'Step 2 Size'
                        T = table([meanStep2; stdStep2]);
                        colNames = {'step 2 sizes (nm)','probability'...
                            'fit sizes (nm)','fit probability'};
                        rowNames = {'mean (nm)';'standard deviation (nm)'};
                    case 'Total Step Size'
                        T = table([meanTotalStep; stdTotalStep]);
                        colNames = {'total step sizes (nm)','probability'...
                            'fit sizes (nm)','fit probability'};
                        rowNames = {'mean (nm)';'standard deviation (nm)'};
                    case 'Event Duration / Force on Selected Bead'
                        T = table; % no parameters to output
                        colNames = {'force on selected bead (pN)','log(event duration)'};
                    otherwise
                        T = table(kDur);
                        colNames = {'event durations (s)','probability'...
                            'fit durations (s)','fit probability'};
                        rowNames = {'k'};
                end

                % Write the data to the specified .xlsx file.
                writetable(table(colNames),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(1) '1']);
                writetable(table(miscData),...
                    [path name '.xlsx'],...
                    'WriteVariableNames', 0,...
                    'Range', [app.getExcelCol(1) '2']);
                if ~isempty(miscFit)
                    writetable(table(miscFit),...
                        [path name '.xlsx'],...
                        'WriteVariableNames', 0,...
                        'Range', [app.getExcelCol(3) '2']);
                end
                if ~isempty(T)
                    writetable(table(rowNames),...
                        [path name '.xlsx'],...
                        'WriteVariableNames', 0,...
                        'Range', [app.getExcelCol(6) '1']);
                    writetable(T,...
                        [path name '.xlsx'],...
                        'WriteVariableNames', 0,...
                        'Range', [app.getExcelCol(7) '1']);
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME, sprintf('An unexpected error occurred while saving:\n\n%s', ME.message));
                app.restore(store);
            end
        end
        
        function chooseDataCallback(app, src, ~, ax)
            % Called when user clicks 'Choose Data' button on the UIContextMenu of ax.
            % Opens a dialog box allowing them to switch which data is shown in the axes.
            
            switch ax.UserData
                case 'miscAxes'
                    backup = app.storeBackupTab1;
                case 'combinedMiscAxes'
                    backup = app.storeBackupTab2;
            end
            
            try
                if app.misc.disabled
                    return
                end

                store = app.disable;

                % Determine which function to call.
                switch ax.UserData
                    case 'miscAxes'
                        % Only continue if selected button was not already selected.
                        if strcmp(app.analyze.misc.data, src.Text)
                            app.restore(store);
                            return;
                        end
                        app.analyze.misc.data = src.Text;
                        app.plotMisc;
                    case 'combinedMiscAxes'
                        % Only continue if selected button was not already selected.
                        if strcmp(app.analyze2.misc.data, src.Text)
                            app.restore(store);
                            return;
                        end
                        app.analyze2.misc.data = src.Text;
                        app.plotMisc2;
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                switch ax.UserData
                    case 'miscAxes'
                        app.restoreBackupTab1(backup);
                    case 'combinedMiscAxes'
                        app.restoreBackupTab2(backup);
                end
            end
        end
        
        function updateCurrentTask(app, s)
            % Updates the currentTask text to display s, and repositions the text based on
            % the size of s.
            
            app.load.currentTask.String = s;
            width = min([app.load.taskMaxWidth, app.load.currentTask.Extent(3)]);
            app.load.currentTask.Position(1) = 0.995 - width;
            app.load.currentTask.Position(3) = width;
            
            app.load2.currentTask.String = s;
            width = min([app.load2.taskMaxWidth, app.load2.currentTask.Extent(3)]);
            app.load2.currentTask.Position(1) = 0.995 - width;
            app.load2.currentTask.Position(3) = width;
            
            drawnow;
        end
        
        function enableable = getEnableable(app)
            % Returns all 'enableable' uicontrols which should be disabled while the
            % program is busy.
            
            enableable = [app.load.inputBtn,...
                app.analyze.individual.showRemovedBtn,...
                app.analyze.individual.currentEventTxt,...
                app.analyze.individual.currentEventEdit,...
                app.analyze.individual.ofNumTxt,...
                app.analyze.individual.prevEventBtn,...
                app.analyze.individual.nextEventBtn,...
                app.analyze.individual.removeBtn,...
                app.analyze.individual.pickChngptBtn,...
                app.analyze.individual.resetBtn,...
                app.analyze.individual.snap,...
                app.controls.deselect.deselectBtn,...
                app.controls.generateCovBtn,...
                app.controls.findEventsBtn,...
                app.controls.useChangepointBtn,...
                app.controls.toggleEventsBtn,...
                app.controls.flipDataBtn,...
                app.controls.peakToPeakBtn,...
                app.controls.thresholdBtn,...
                app.controls.automaticBtn,...
                app.controls.manualBtn,...
                app.controls.setManualBtn,...
                app.controls.peak1Edit,...
                app.controls.peak2Edit,...
                app.controls.minEdit,...
                app.controls.beadABtn,...
                app.controls.beadBBtn,...
                app.controls.bothBeadsBtn,...
                app.controls.windows.covwindowEdit,...
                app.controls.windows.covsmoothEdit,...
                app.controls.windows.interactBtn,...
                app.controls.windows.resetBtn,...
                app.controls.cutoffs.minDurEdit,...
                app.controls.cutoffs.minSepEdit,...
                app.controls.cutoffs.resetBtn,...
                app.excel.listbox,...
                app.excel.createBtn,...
                app.excel.findBtn,...
                app.excel.addBtn,...
                app.controls.windows.covwindowTxt,...
                app.controls.windows.covsmoothTxt,...
                app.controls.cutoffs.minDurTxt,...
                app.controls.cutoffs.minSepTxt,...
                app.controls.peak1Txt,...
                app.controls.peak2Txt,...
                app.controls.minTxt,...
                app.load2.inputBtn,...
                app.load2.table,...
                app.load2.beadABtn,...
                app.load2.beadBBtn,...
                app.load2.bothBeadsBtn];
        end
        
        function store = disable(app, varargin)
            % Disables all 'enableable' uicontrols and returns an array of uicontrols
            % which were previously enabled.
            
            app.misc.disabled = true;
            enableable = app.getEnableable;
            store = zeros(1, length(enableable));
            for i = 1:length(enableable)
                store(i) = strcmp(enableable(i).Enable, 'on');
                enableable(i).Enable = 'off';
            end
            if nargin == 2 && strcmp(varargin{1}, '-toolbars')
                % If '-toolbars' is specified, return early.
                return
            else
                % Otherwise, disable the toolbars.
                if app.misc.axToolbars
                    % R2018b or later (axes toolbars):
                    app.setAxToolbar(app.load.main.mainAxes, 'off');
                    app.setAxToolbar(app.load.cov.covAxes, 'off');
                    app.setAxToolbar(app.analyze.covHist.covHistAxes, 'off');
                    app.setAxToolbar(app.analyze.misc.miscAxes, 'off');
                    app.setAxToolbar(app.analyze.ensemble.ensembleAxes, 'off');
                    app.setAxToolbar(app.analyze.individual.likelihoodAxes, 'off');
                    app.setAxToolbar(app.analyze.individual.eventAxes, 'off');
                    app.setAxToolbar(app.analyze2.misc.miscAxes, 'off');
                    app.setAxToolbar(app.analyze2.ensemble.ensembleAxes, 'off');
                else
                    % pre-R018b (menu items):
                    set(app.misc.menus, 'Enable', 'off')
                end
            end
        end
        
        function addToNew(app, on, objs)
            % Updates app.misc.toBeEnabled to reflect whether uicontrols contained in objs
            % should be enabled (on = true) or disabled (on = false) after program is
            % restored.
            
            if on
                toBeEnabled = 1;
            else
                toBeEnabled = -1;
            end
            for i = 1:length(objs)
                app.misc.toBeEnabled(objs(i) == app.getEnableable) = toBeEnabled;
            end
        end
        
        function restore(app, store)
            % Resets all uicontrols to their original state prior to disable(), as
            % determined by values in store, or updates their state based on
            % app.misc.toBeEnabled.
            
            enableable = app.getEnableable;
            for i = 1:length(enableable)
                if app.misc.toBeEnabled(i) == 1
                    % uicontrol should be newly enabled.
                    enableable(i).Enable = 'on';
                elseif app.misc.toBeEnabled(i) == -1
                    % uicontrol should be newly disabled.
                    enableable(i).Enable = 'off';
                elseif store(i)
                    % uicontrol should enable if it was enabled previously.
                    enableable(i).Enable = 'on';
                end
            end
            if app.misc.axToolbars
                app.setAxToolbar(app.load.main.mainAxes, 'on');
                app.setAxToolbar(app.load.cov.covAxes, 'on');
                app.setAxToolbar(app.analyze.covHist.covHistAxes, 'on');
                app.setAxToolbar(app.analyze.misc.miscAxes, 'on');
                app.setAxToolbar(app.analyze.ensemble.ensembleAxes, 'on');
                app.setAxToolbar(app.analyze.individual.likelihoodAxes, 'on');
                app.setAxToolbar(app.analyze.individual.eventAxes, 'on');
                app.setAxToolbar(app.analyze2.misc.miscAxes, 'on');
                app.setAxToolbar(app.analyze2.ensemble.ensembleAxes, 'on');
            else
                set(app.misc.menus, 'Enable', 'on');
            end
            
            % Reset current task.
            app.updateCurrentTask('');
            
            % Reset app.misc.toBeEnabled.
            app.misc.toBeEnabled = zeros(length(app.getEnableable), 1);
            
            app.misc.disabled = false;
        end
        
        function setupAxToolbar(app, ax)
            % Sets up ax toolbar to have only 'export', 'pan', 'zoom', and 'restore'
            % buttons and establishes SelectionChangedFcn.
            
            if app.misc.axToolbars
                axtoolbar(ax, {'export','pan','zoomin','zoomout','rotate','restoreview'},...
                    'SelectionChangedFcn', @app.axToolbarSelectionChangedCallback,...
                    'Visible', 'off');
                try
                    disableDefaultInteractivity(ax);
                catch
                end
            end
        end
        
        function axToolbarSelectionChangedCallback(app, ~, event, store)
            % Detects if user has clicked an axis toolbar button and disables app if so.
            % This is only called for state buttons, i.e. zoom, pan, rotate.
            
            try
                state = event.Selection.Value;
                if strcmp(state, 'on')
                    % If the selected button is on, update the current task message.
                    app.updateCurrentTask(sprintf('Editing plots in %s mode.', event.Selection.Tooltip));
                    if nargin < 4
                        % If no store has already been passed to this callback (i.e. app isn't
                        % already disabled), disable the app and pass store to this callback.
                        store = app.disable('-toolbars');
                        updateSelectionChangedFcn(app.load.main.mainAxes, true);
                        updateSelectionChangedFcn(app.load.cov.covAxes, true);
                        updateSelectionChangedFcn(app.analyze.covHist.covHistAxes, true);
                        updateSelectionChangedFcn(app.analyze.misc.miscAxes, true);
                        updateSelectionChangedFcn(app.analyze.ensemble.ensembleAxes, true);
                        updateSelectionChangedFcn(app.analyze.individual.likelihoodAxes, true);
                        updateSelectionChangedFcn(app.analyze.individual.eventAxes, true);
                        updateSelectionChangedFcn(app.analyze2.misc.miscAxes, true);
                        updateSelectionChangedFcn(app.analyze2.ensemble.ensembleAxes, true);
                    end
                elseif strcmp(state, 'off')
                    % If the selected button is off, reset the callback and restore the app.
                    updateSelectionChangedFcn(app.load.main.mainAxes, false);
                    updateSelectionChangedFcn(app.load.cov.covAxes, false);
                    updateSelectionChangedFcn(app.analyze.covHist.covHistAxes, false);
                    updateSelectionChangedFcn(app.analyze.misc.miscAxes, false);
                    updateSelectionChangedFcn(app.analyze.ensemble.ensembleAxes, false);
                    updateSelectionChangedFcn(app.analyze.individual.likelihoodAxes, false);
                    updateSelectionChangedFcn(app.analyze.individual.eventAxes, false);
                    updateSelectionChangedFcn(app.analyze2.misc.miscAxes, false);
                    updateSelectionChangedFcn(app.analyze2.ensemble.ensembleAxes, false);
                    app.restore(store);
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                try
                    app.restore(store);
                catch
                end
            end
            
            function updateSelectionChangedFcn(ax, on)
                if on
                    ax.Toolbar.SelectionChangedFcn = {@app.axToolbarSelectionChangedCallback, store};
                else
                    ax.Toolbar.SelectionChangedFcn = @app.axToolbarSelectionChangedCallback;
                end
            end
        end
        
        function setAxToolbar(app, ax, on)
            % Sets ax toolbar visibility to value of on, given that ax has a toolbar and
            % is itself visible.
            
            if app.misc.axToolbars
                if strcmp(ax.Visible, 'off')
                    ax.Toolbar.Visible = 'off';
                else
                    ax.Toolbar.Visible = on;
                end
            end
        end
        
        function menuSelectedFcn(app, src, ~, store)
            % Called whenever pre-R2018b user selects any menu item in the menu bar. Calls
            % the default callback function for that menu item and disables or restores
            % the app if needed.
            
            try
                % Call the default callback method for src.
                switch src.Tag
                    case 'figMenuPan'
                        toolsmenufcn(gcbf, 'Pan');
                    case 'figMenuRotate3D'
                        toolsmenufcn(gcbf, 'Rotate');
                    case 'figMenuZoomIn'
                        toolsmenufcn(gcbf, 'ZoomIn');
                    case 'figMenuZoomOut'
                        toolsmenufcn(gcbf, 'ZoomOut');
                    case 'figMenuResetView'
                        toolsmenufcn(gcbf, 'ResetView');
                    case 'figMenuFileSaveAs'
                        filemenufcn(gcbf, 'FileSaveAs');
                end

                idx = cellfun(@(c) strcmp(c, 'on'), {app.misc.submenus.Checked});
                if any(idx)
                    % If any of the menu items which can be checked are checked, update the
                    % current task message.
                    app.updateCurrentTask(sprintf('Editing plots in %s mode.', strrep(app.misc.submenus(idx).Text,'&','')));
                    if nargin < 4
                        % If no store has already been passed to this callback, then disable
                        % the app and pass store to this callback for each menu item.
                        store = app.disable('-toolbars');
                        set(app.misc.submenus, 'MenuSelectedFcn', {@app.menuSelectedFcn store});
                    end
                elseif nargin == 4
                    % If no menu items are checked and store has been passed to this callback,
                    % restore the app and reset each menu item's callback.
                    set(app.misc.submenus, 'MenuSelectedFcn', @app.menuSelectedFcn);
                    app.restore(store);
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                try
                    app.restore(store);
                catch
                end
            end
        end
        
        function skipBtnCallback(app, ~, ~)
            % Called when user clicks skip button when program is calculating ensemble
            % averages.
            
            app.misc.skip = true;
            app.load.skipBtn.Enable = 'off';
            app.load2.skipBtn.Enable = 'off';
        end
        
        function postSetCallback(app, src, ~)
            % Called just after the value of src is set.
            
            app.misc.listenerError = [];
            try
                switch src.Name
                    case 'showEvents'
                        % showEvents has changed, so event patches need to be updated.
                        app.updateEventVisNColor;
                    case 'showRemoved'
                        % showRemoved has changed, so current event, individual
                        % event/likelihood plots, and event patches possibly need to be
                        % updated.
                        if ~isempty(app.activeEvents)
                            % If events have been found...
                            ce = app.analyze.individual.storeCurrentEvent; % Previous current event.
                            switch app.showRemoved
                                case 0
                                    % Removed events should now be hidden, so in case the
                                    % previous current event is no longer valid, the new
                                    % current event should be its nearest valid neighbor.
                                    k = app.closestTrue(app.activeEvents);
                                    j = sum(app.activeEvents(1:k(ce))); % New current event.
                                    N = sum(app.activeEvents);
                                case 1
                                    % Removed events were invisible and are now visible, so
                                    % previous current event must have been valid and is
                                    % definitely still valid.
                                    j = ce; % New current event.
                                    N = length(app.activeEvents);
                            end
                            app.changeCurrentEvent(j); % Update current event.
                            if ~app.activeEvents(ce)
                                % If the previous current event was inactive, plots might need
                                % to be updated.
                                app.plotIndividualEvent;
                                app.updateEventVisNColor; % Update event patches.
                            end
                            app.analyze.individual.ofNumTxt.String = ['of ', num2str(N)];
                        end
                        app.analyze.individual.showRemovedBtn.Value = app.showRemoved;
                    case 'activeEvents'
                        % activeEvents has changed, so current event, individual
                        % event/likelihood plots, event patches, misc axes, and ensemble
                        % average plot possibly need to be updated.
                        if ~isempty(app.activeEvents)
                            % If there are any events at all...
                            if ~any(app.activeEvents)
                                % There are events, but they have all been removed, so set
                                % showRemoved to true, hide misc axes and ensemble axes, and
                                % tell user what is going on.
                                app.analyze.individual.storeCurrentEvent = 1;
                                app.showRemoved = true;
                                if ~isempty(app.misc.listenerError)
                                    rethrow(app.misc.listenerError)
                                end
                                app.noActiveEventsReset;
                                warndlg('All of the potential events have been removed.')
                            else
                                % There are events which have not been removed, so update
                                % plots.
                                app.updateEventVisNColor; % First pass to hide removed events.
                                app.plotMisc;
                                [app.analyze.ensemble.start, app.analyze.ensemble.end] = app.calcEnsemble;
                                app.plotEnsemble;
                                if app.misc.optim
                                    [app.analyze.ensemble.startFit, app.analyze.ensemble.endFit,...
                                        app.analyze.ensemble.p] = app.calcEnsembleFits;
                                    app.plotEnsembleFits;
                                end
                            end
                            app.addToNew(any(app.activeEvents), app.analyze.individual.showRemovedBtn); % Only let user change showRemoved if there are active events.
                            ce = app.analyze.individual.storeCurrentEvent; % Previous current event.
                            switch app.showRemoved
                                case 0
                                    if isempty(ce)
                                        % If there were previously no events, then just show
                                        % the first event.
                                        app.changeCurrentEvent(1);
                                        app.plotIndividualEvent;
                                    else
                                        k = app.closestTrue(app.activeEvents);
                                        j = sum(app.activeEvents(1:k(ce))); % New current event.
                                        app.changeCurrentEvent(j);
                                        if ~app.activeEvents(ce)
                                            % If the previous current event was inactive,
                                            % plots might need to be updated.
                                            app.plotIndividualEvent;
                                        end
                                    end
                                    N = sum(app.activeEvents);
                                case 1
                                    N = length(app.activeEvents);
                            end
                            app.analyze.individual.totalPotentialTxt.String = ['Number of Potential Events: ', num2str(length(app.activeEvents))];
                            app.analyze.individual.totalActiveTxt.String = ['Number of Analyzed Events: ', num2str(sum(app.activeEvents))];
                            app.analyze.individual.ofNumTxt.String = ['of ', num2str(N)];
                            app.updateCheckboxAndLabel; % Update removed checkbox and label.
                            app.updateEventVisNColor; % Second pass to show the current event.
                        end
                    case 'allEvents'
                        % allEvents has changed, so events should be plotted and activeEvents
                        % should be updated.
                        if ~isempty(app.allEvents)
                            % If there are any events at all...
                            app.plotEvents; % ...update the event plots...
                            app.analyze.individual.storeCurrentEvent = app.getCurrentEventAsIndex;
                            app.activeEvents = [];
                            app.pruneEvents; % ...and update activeEvents.
                        else
                            app.analyze.individual.removedByUser = [];
                            app.analyze.individual.includedByUser = [];
                        end
                    case 'flipped'
                        % flipped has changed, so main plot, ensemble average plot, and
                        % individual event plot need to be updated.
                        if strcmp(app.load.main.mainAxes.Visible, 'on')
                            % Only update if showing.
                            app.plotGlobal;
                        end
                        if strcmp(app.analyze.ensemble.ensembleAxes.Visible, 'on')
                            % Only update if showing.
                            app.plotEnsemble;
                            app.plotEnsembleFits;
                        end
                        if strcmp(app.analyze.individual.eventAxes.Visible, 'on')
                            % Only update if showing.
                            app.initializeIndividualEvent;
                            app.plotIndividualEvent;
                        end
                    case {'minSep', 'minDur'}
                        % minSep or minDur have changed, so activeEvents might need to be
                        % updated.
                        if ~isempty(app.allEvents)
                            % If there are any events at all...
                            app.analyze.individual.storeCurrentEvent = app.getCurrentEventAsIndex;
                            app.pruneEvents; % ...update activeEvents.
                        end
                    case {'autop1', 'autop2', 'automin', 'manp1', 'manp2', 'manmin', 'usemin', 'useman'}
                        % Any of the peaks/min, usemin, or useman have changed, so covariance
                        % and covariance histogram plots might need to be updated and events
                        % might be outdated.

                        % Obtain the correct peak/min values based on manual vs automatic.
                        if app.useman
                            m = app.manmin;
                            p1 = app.manp1;
                            p2 = app.manp2;
                        else
                            m = app.automin;
                            p1 = app.autop1;
                            p2 = app.autop2;
                        end

                        % Update edit texts regardless of peak-to-peak vs threshold.
                        app.controls.peak1Edit.String = num2str(p1);
                        app.controls.peak2Edit.String = num2str(p2);
                        app.controls.minEdit.String = num2str(m);

                        % Now take peak-to-peak vs threshold into account.
                        if app.usemin
                            p1 = [];
                            p2 = [];
                        else
                            m = [];
                        end

                        % Update covariance and covariance histogram plots.
                        if strcmp(app.analyze.covHist.covHistAxes.Visible, 'on')
                            app.updateCurrentTask('Plotting the peaks and minimum.');
                        end
                        app.updatePeaksMinPlots(p1, p2, m);
                        drawnow;

                        % As covariance/peaks/min may be different, check if events need to be
                        % updated.
                        app.checkOutdatedEvents;
                        
                        % Update radio buttons.
                        app.controls.automaticBtn.Value = ~app.useman;
                        app.controls.manualBtn.Value = app.useman;
                        app.controls.peakToPeakBtn.Value = ~app.usemin;
                        app.controls.thresholdBtn.Value = app.usemin;
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                % Exceptions thrown within this function are typically caught in some file
                % written by TMW and turned into warnings. Instead, we want to pass the
                % exception back to whichever line of code invoked this function, display
                % it to the user, and terminate whichever function contains that line of
                % code. So, we store it here and then check if app.misc.listenerError is
                % empty whenever we change an obversable property.
                app.misc.listenerError = ME;
            end
        end
        
        function varargout = errordlgTrace(app, msg, trace)
            % Displays msg in a dialog box which mimics the default error dialog box. Lets
            % user click a button to show trace in an edit text. Created with the help of
            % errordlg and questdlg.
            
            % Get texts.
            msg = cellstr(msg);
            trace = cellstr(trace);
            title = getString(message('MATLAB:uistring:popupdialogs:ErrorDialogTitle'));
            btn1 = getString(message('MATLAB:uistring:popupdialogs:OK'));
            btn2 = 'See Trace';

            % Set up initial widths, heights, offsets.
            defOffset = 7;
            iconWidth = 32 * 72/get(groot, 'ScreenPixelsPerInch');
            iconHeight = 32 * 72/get(groot, 'ScreenPixelsPerInch');
            iconXOffset = defOffset;

            figPos = get(0, 'DefaultFigurePosition');
            figWidth = 190;
            figHeight = 50;
            figPos(3:4) = [figWidth, figHeight]; % May change.

            btnWidth = 50;
            btnHeight = 17;
            btnYOffset = defOffset;

            msgTxtXOffset = iconXOffset+iconWidth+defOffset;
            msgTxtYOffset = defOffset+btnYOffset+btnHeight; % May change.
            msgTxtWidth = figWidth-2*defOffset-iconWidth; % May change.
            msgTxtHeight = figHeight-defOffset-msgTxtYOffset; % May change.

            % Create main window.
            fig = dialog('Visible', 'off',...
                'Name', title,...
                'Units', 'points',...
                'Position', figPos);

            % Add buttons - x offsets pending as figure size not yet established.
            font.FontUnits = 'points';
            font.FontSize = get(0, 'FactoryUicontrolFontSize');
            font.FontName = get(0, 'FactoryUicontrolFontName');
            font.FontWeight = get(fig, 'DefaultUicontrolFontWeight');
            btnHandle{1} = uicontrol(fig,...
                font,...
                'Style', 'pushbutton',...
                'Units', 'points',...
                'Position', [0, btnYOffset, btnWidth, btnHeight],...
                'KeyPressFcn', @doControlKeyPress,...
                'Callback', @(~,~) delete(gcbf),...
                'String', btn1,...
                'HorizontalAlignment', 'center');
            btnHandle{2} = uicontrol(fig,...
                font,...
                'Style', 'pushbutton',...
                'Units', 'points',...
                'Position', [0, btnYOffset, btnWidth, btnHeight],...
                'KeyPressFcn', @doControlKeyPress,...
                'Callback', @seeFullErrorMessage,...
                'String', btn2,...
                'HorizontalAlignment', 'center');

            % Wrap message, with max 75 characters per line.
            stFont = font;
            stFont.FontWeight = get(fig, 'DefaultTextFontWeight');
            msgHandle = uicontrol(fig,...
                stFont,...
                'Style', 'text',...
                'Units', 'points',...
                'Position', [0, 0, 0, 0],...
                'String', {' '},...
                'HorizontalAlignment', 'left',...
                'BackgroundColor', fig.Color,...
                'ForegroundColor', 'k');
            [wrapString, newMsgTxtPos] = textwrap(msgHandle, msg, 75);

            % Wrap trace similarly.
            [wrapStringTrace, newMsgTxtPosTrace] = textwrap(msgHandle, trace, 75);
            delete(msgHandle);

            % Display message and determine new width and height.
            msgAxes = axes('Parent', fig, 'Position', [0 0 1 1], 'Visible', 'off');
            textHandle = text('Parent', msgAxes,...
                'Units', 'points',...
                'String', wrapString,...
                'Color', btnHandle{1}.ForegroundColor,...
                stFont,...
                'HorizontalAlignment', 'left',...
                'VerticalAlignment', 'bottom',...
                'Interpreter', 'none');
            msgTxtWidth = max([msgTxtWidth, newMsgTxtPos(3), textHandle.Extent(3)]);
            msgTxtHeight = max([msgTxtHeight, newMsgTxtPos(4), textHandle.Extent(4)]);
            
            % Similar for trace, using an edit text so user can copy.
            traceBox = uicontrol(fig,...
                stFont,...
                'Style', 'edit',...
                'Units', 'points',...
                'Max', 2,...
                'String', wrapStringTrace,...
                'Visible', 'off',...
                'HorizontalAlignment', 'left');
            msgTxtWidthTrace = max([msgTxtWidth, newMsgTxtPosTrace(3), traceBox.Extent(3)]);
            msgTxtHeightTrace = min([max([msgTxtHeight, newMsgTxtPosTrace(4), traceBox.Extent(4)]), 350]);
            msgTxtWidthTrace = msgTxtWidthTrace + ceil(20 * 72/get(groot, 'ScreenPixelsPerInch')); % Account for scrollbar.
            
            textHandle.String = wrapString;
            
            % Vertically center icon and message.
            if iconHeight > msgTxtHeight
                iconYOffset = btnYOffset+btnHeight+defOffset;
                msgTxtYOffset = iconYOffset+(iconHeight-msgTxtHeight)/2;
                figPos(4) = iconYOffset+iconHeight+defOffset;
            else
                msgTxtYOffset = btnYOffset+btnHeight+defOffset;
                iconYOffset = msgTxtYOffset+(msgTxtHeight-iconHeight)/2;
                figPos(4) = msgTxtYOffset+msgTxtHeight+defOffset;
            end

            % Main window width is determined either by buttons or by message.
            figPos(3) = max(2*(btnWidth+defOffset)+defOffset, msgTxtXOffset+msgTxtWidth+defOffset);
            set(fig, 'Position', reposition(figPos));

            % Button x offsets can now be set as figure size is final.
            btnXOffset = [(figPos(3)-defOffset)/2-btnWidth, (figPos(3)+defOffset)/2];
            btnHandle{1}.Position(1) = btnXOffset(1);
            btnHandle{2}.Position(1) = btnXOffset(2);

            % Place message.
            textHandle.Position = [msgTxtXOffset, msgTxtYOffset, 0];

            % Plot error icon.
            iconAxes = axes('Parent', fig,...
                'Units', 'points',...
                'Position', [iconXOffset, iconYOffset, iconWidth, iconHeight]);
            [iconData, alphaData] = matlab.ui.internal.dialog.DialogUtils.imreadDefaultIcon('error');
            img = image('CData', iconData, 'AlphaData', alphaData, 'Parent', iconAxes);
            set(iconAxes,...
                'Visible', 'off',...
                'YDir', 'reverse',...
                'XLim', img.XData+[-0.5 0.5], ...
                'YLim', img.YData+[-0.5 0.5]);
            
            % Take 2, now with trace message:
            
            % Vertically center icon and message.
            if iconHeight > msgTxtHeightTrace
                iconYOffset = btnYOffset+btnHeight+defOffset;
                msgTxtYOffset = iconYOffset+(iconHeight-msgTxtHeightTrace)/2;
                figPos(4) = iconYOffset+iconHeight+defOffset;
            else
                msgTxtYOffset = btnYOffset+btnHeight+defOffset;
                iconYOffset = msgTxtYOffset+(msgTxtHeightTrace-iconHeight)/2;
                figPos(4) = msgTxtYOffset+msgTxtHeightTrace+defOffset;
            end
            
            % Main window width is determined either by buttons or by message.
            figPos(3) = max(2*(btnWidth+defOffset)+defOffset, msgTxtXOffset+msgTxtWidthTrace+defOffset);
            
            % Button x offsets can now be set as figure size is final.
            btnXOffset = [(figPos(3)-defOffset)/2-btnWidth, (figPos(3)+defOffset)/2];
            
            % Place edit text.
            traceBox.Position = [msgTxtXOffset, msgTxtYOffset, msgTxtWidthTrace, msgTxtHeightTrace];

            % Make sure window is on screen, make visible, give focus to OK button.
            movegui(fig)
            fig.Visible = 'on';
            uicontrol(btnHandle{1});
            drawnow

            % Return handle to figure if requested, so you can call uitwait() on this.
            if nargout == 1
                varargout{1} = fig;
            end

            function seeFullErrorMessage(varargin)
                % Reposition figure.
                fig.Position(3:4) = figPos(3:4);
                movegui(fig);
                
                % Reposition buttons.
                btnHandle{1}.Position(1) = btnXOffset(1);
                btnHandle{2}.Position(1) = btnXOffset(2);
                
                % Show stack trace instead of basic error.
                textHandle.Visible = 'off';
                traceBox.Visible = 'on';
                traceBox.String = wrapStringTrace;
                
                % Reposition icon.
                iconAxes.Position = [iconXOffset, iconYOffset, iconWidth, iconHeight];
            end
            
            function doControlKeyPress(~, evd)
                switch evd.Key
                    case 'return'
                        evd.Source.Callback();
                    case 'escape'
                        delete(fig)
                end
            end
            
            function size = reposition(size)
                ogunits = app.misc.fig.Units;
                app.misc.fig.Units = 'points';
                size(1) = app.misc.fig.Position(1) + 1/2*(app.misc.fig.Position(3) - size(3));
                size(2) = app.misc.fig.Position(2) + 2/3*(app.misc.fig.Position(4) - size(4));
                app.misc.fig.Units = ogunits;
            end
        end
        
        function unexpectedError(app, ME, str)
            % Displays the exception contained in ME with a modal error dialog box.
            
            if nargin < 3
                % Default message.
                str = sprintf('An unexpected error occurred:\n\n%s', ME.message);
            end
            try
                % First, try to use custom error dialog box.
                uiwait(app.errordlgTrace(str, getReport(ME, 'extended', 'hyperlinks', 'off')));
            catch
                % If unable to do so, for whatever reason, use TMW's error dialog box.
                uiwait(errordlg(str, 'Error Dialog', 'modal'));
            end
        end
    end
    
    % Private methods for simulated input
    methods (Access = private)
        
        function events = getRealEvents(app)
            % Returns the real events within the simulated data.
            
            events = app.misc.realEvents;
        end
        
        function dur = getRealDur(app)
            % Returns the durations of the real events.
            
            events = app.getRealEvents;
            if ~isempty(events)
                dur = events(:,2) - events(:,1);
            else
                dur = [];
            end
        end
        
        function sep = getRealSep(app)
            % Returns the number of points between consecutive events for the real events.
            
            events = app.getRealEvents;
            if ~isempty(events)
                sep = [events(:,1); length(app.getTrimBeads)] - [1; events(:,2)];
            else
                sep = [];
            end
        end
        
        function [pos, steps] = getRealPosNStepSizes(app, whichBead)
            % Returns the average position of a bead just before and just after binding and
            % unbinding, as well as the estimated step sizes, for the real events.
            
            events = app.getRealEvents;
            [A, B] = app.getActiveBeads;
            if app.flipped
                A = -A;
                B = -B;
            end
            if strcmp(whichBead, 'B')
                bead = B;
            else
                bead = A;
            end
            N = size(events, 1);
            M = numel(bead);
            posBefB = zeros(N, 1);
            posAftB = zeros(N, 1);
            posBefU = zeros(N, 1);
            posAftU = zeros(N, 1);
            events = [0 0; events; M, M];
            for i = 2:N+1
                posBefB(i-1) = mean(bead(max([events(i-1,2)+1, events(i,1)-app.data.fs/100]):events(i,1)));
                posAftB(i-1) = mean(bead(events(i,1)+1:min([events(i,1)+1+app.data.fs/100, events(i,2)])));
                posBefU(i-1) = mean(bead(max([events(i,2)-app.data.fs/100, events(i,1)+1]):events(i,2)));
                posAftU(i-1) = mean(bead(events(i,2)+1:min([events(i,2)+1+app.data.fs/100, events(i+1,1)])));
            end
            step1 = posAftB - posBefB;
            totalstep = posBefU - posAftU;
            step2 = totalstep - step1;
            pos = [posBefB...
                posAftB...
                posBefU...
                posAftU];
            steps = [step1...
                step2...
                totalstep];
        end
        
        function plotRealEvents(app, events)
            % Plots real events on the main axes and covariance axes in load panel of tab 1.
            
            app.updateCurrentTask('Plotting the real events.');
            
            if ~app.misc.simInput
                return
            end
            
            % Hide all previous real event patches.
            delete(app.misc.simPatches);
            app.misc.simPatches = gobjects(0);
            
            % Plot new real event patches.
            limyC = ylim(app.load.cov.covAxes);
            limyM = ylim(app.load.main.mainAxes);
            liminitC = limyC; liminitM = limyM;
            limyC(1) = limyC(1) - 0.1*(limyC(2) - limyC(1)); % Each patch will span from 10% below...
            limyC(2) = limyC(2) + 0.1*(limyC(2) - limyC(1)); % ...to 10% above current ylim.
            limyM(1) = limyM(1) - 0.1*(limyM(2) - limyM(1)); % Repeat for main axes.
            limyM(2) = limyM(2) + 0.1*(limyM(2) - limyM(1));
            fileTime = app.data.time;
            for nevent = 1:size(events, 1)
                % For each event...
                event = events(nevent,:);
                app.misc.simPatches(end+1) = patch(app.load.cov.covAxes,...
                    'XData', [fileTime(event(1)+1) fileTime(event(2)) fileTime(event(2)) fileTime(event(1)+1)],...
                    'YData', [limyC(1) limyC(1) limyC(2) limyC(2)],...
                    'FaceColor', app.misc.colors.simPatches,...
                    'FaceAlpha', 0.3,...
                    'EdgeColor', 'none',...
                    'HitTest', 'off'); % ...plot a corresponding patch.
                app.misc.simPatches(end+1) = patch(app.load.main.mainAxes,...
                    'XData', [fileTime(event(1)+1) fileTime(event(2)) fileTime(event(2)) fileTime(event(1)+1)],...
                    'YData', [limyM(1) limyM(1) limyM(2) limyM(2)],...
                    'FaceColor', app.misc.colors.simPatches,...
                    'FaceAlpha', 0.3,...
                    'EdgeColor', 'none',...
                    'HitTest', 'off'); % ...plot a corresponding patch.
            end
            ylim(app.load.cov.covAxes, liminitC);
            ylim(app.load.main.mainAxes, liminitM);
            
            drawnow;
        end
    end
    
    % Private tab 1 methods
    methods (Access = private)
        
        function backup = storeBackupTab1(app)
            % Copies by value - not by reference - all tab 1 properties which could
            % potentially change.
            
            backup.misc.fig.Position = app.misc.fig.Position;
            backup.misc.fig.WindowButtonMotionFcn = app.misc.fig.WindowButtonMotionFcn;
            backup.misc.fig.WindowButtonDownFcn = app.misc.fig.WindowButtonDownFcn;
            backup.misc.fig.WindowKeyPressFcn = app.misc.fig.WindowKeyPressFcn;
            backup.misc.fig.SizeChangedFcn = app.misc.fig.SizeChangedFcn;
            backup.misc.storeClickCallback = app.misc.storeClickCallback;
            backup.misc.storeMotionCallback = app.misc.storeMotionCallback;
            backup.misc.fig.Pointer = app.misc.fig.Pointer;
            
            backup.load.inputBtn = storeButton(app.load.inputBtn);
            backup.load.path = app.load.path;
            backup.load.name = app.load.name;
            backup.load.ext = app.load.ext;
            backup.load.header = app.load.header;
            backup.data.origData = app.data.origData;
            backup.data.origBeads = app.data.origBeads;
            backup.data.time = app.data.time;
            backup.load.CALA = app.load.CALA;
            backup.load.CALATxt = storeText(app.load.CALATxt);
            backup.load.CALB = app.load.CALB;
            backup.load.CALBTxt = storeText(app.load.CALBTxt);
            backup.load.KA = app.load.KA;
            backup.load.KATxt = storeText(app.load.KATxt);
            backup.load.KB = app.load.KB;
            backup.load.KBTxt = storeText(app.load.KBTxt);
            backup.data.fs = app.data.fs;
            backup.load.fsTxt = storeText(app.load.fsTxt);
            backup.load.filenameTxt = storeText(app.load.filenameTxt);
            backup.load.skipBtn = storeButton(app.load.skipBtn);
            backup.load.skipBtn.Position = app.load.skipBtn.Position;
            backup.load.taskMaxWidth = app.load.taskMaxWidth;
            backup.misc.simInput = app.misc.simInput;
            backup.misc.realEvents = app.misc.realEvents;
            
            backup.load.main.mainAxes = storeAxes(app.load.main.mainAxes);
            backup.load.main.mainAxes.XLabel.Visible = app.load.main.mainAxes.XLabel.Visible;
            backup.load.main.mainAxes.XTickMode = app.load.main.mainAxes.XTickMode;
            backup.load.main.mainAxes.XTickLabelMode = app.load.main.mainAxes.XTickLabelMode;
            backup.load.main.ylim = app.load.main.ylim;
            backup.load.main.bead1Plot = storePlot(app.load.main.bead1Plot);
            backup.load.main.bead2Plot = storePlot(app.load.main.bead2Plot);
            backup.load.main.patches = storePatch(app.load.main.patches);
            backup.misc.simPatches = storePatch(app.misc.simPatches);
            backup.load.main.label = storeText(app.load.main.label);
            backup.load.main.prevIn = app.load.main.prevIn;
            
            backup.load.cov.covAxes = storeAxes(app.load.cov.covAxes);
            backup.load.cov.covPlot = storePlot(app.load.cov.covPlot);
            backup.load.cov.peak1Plot = storePlot(app.load.cov.peak1Plot);
            backup.load.cov.peak2Plot = storePlot(app.load.cov.peak2Plot);
            backup.load.cov.minPlot = storePlot(app.load.cov.minPlot);
            backup.load.cov.patches = storePatch(app.load.cov.patches);
            backup.load.cov.cov = app.load.cov.cov;
            
            backup.analyze.covHist.covHistAxes = storeAxes(app.analyze.covHist.covHistAxes);
            backup.analyze.covHist.covHistPlot.Data = app.analyze.covHist.covHistPlot.Data;
            backup.analyze.covHist.covHistPlot.Visible = app.analyze.covHist.covHistPlot.Visible;
            backup.analyze.covHist.covHistPlot.NumBins = app.analyze.covHist.covHistPlot.NumBins;
            backup.analyze.covHist.peak1Plot = storePlot(app.analyze.covHist.peak1Plot);
            backup.analyze.covHist.peak2Plot = storePlot(app.analyze.covHist.peak2Plot);
            backup.analyze.covHist.minPlot = storePlot(app.analyze.covHist.minPlot);
            
            backup.analyze.misc.miscAxes = storeAxes(app.analyze.misc.miscAxes);
            backup.analyze.misc.miscAxes.YScale = app.analyze.misc.miscAxes.YScale;
            backup.analyze.misc.miscAxes.YLabel.String = app.analyze.misc.miscAxes.YLabel.String;
            backup.analyze.misc.miscAxes.XLabel.String = app.analyze.misc.miscAxes.XLabel.String;
            backup.analyze.misc.miscPlot = storePlot(app.analyze.misc.miscPlot);
            backup.analyze.misc.miscPlot.Marker = app.analyze.misc.miscPlot.Marker;
            backup.analyze.misc.miscPlot.LineStyle = app.analyze.misc.miscPlot.LineStyle;
            backup.analyze.misc.realPlot = storePlot(app.analyze.misc.realPlot);
            backup.analyze.misc.miscFitPlot = storePlot(app.analyze.misc.miscFitPlot);
            backup.analyze.misc.miscLabel = storeText(app.analyze.misc.miscLabel);
            backup.analyze.misc.data = app.analyze.misc.data;
            backup.analyze.misc.kDur = app.analyze.misc.kDur;
            backup.analyze.misc.meanStep1 = app.analyze.misc.meanStep1;
            backup.analyze.misc.stdStep1 = app.analyze.misc.stdStep1;
            backup.analyze.misc.meanStep2 = app.analyze.misc.meanStep2;
            backup.analyze.misc.stdStep2 = app.analyze.misc.stdStep2;
            backup.analyze.misc.meanTotalStep = app.analyze.misc.meanTotalStep;
            backup.analyze.misc.stdTotalStep = app.analyze.misc.stdTotalStep;
            
            backup.analyze.ensemble.ensembleAxes = storeAxes(app.analyze.ensemble.ensembleAxes);
            backup.analyze.ensemble.ensembleAxes.YLabel.String = app.analyze.ensemble.ensembleAxes.YLabel.String;
            backup.analyze.ensemble.globalFit = app.analyze.ensemble.globalFit;
            backup.analyze.ensemble.start = app.analyze.ensemble.start;
            backup.analyze.ensemble.end = app.analyze.ensemble.end;
            backup.analyze.ensemble.startFit = app.analyze.ensemble.startFit;
            backup.analyze.ensemble.endFit = app.analyze.ensemble.endFit;
            backup.analyze.ensemble.p = app.analyze.ensemble.p;
            backup.analyze.ensemble.startPlot = storePlot(app.analyze.ensemble.startPlot);
            backup.analyze.ensemble.endPlot = storePlot(app.analyze.ensemble.endPlot);
            backup.analyze.ensemble.startFitPlot = storePlot(app.analyze.ensemble.startFitPlot);
            backup.analyze.ensemble.endFitPlot = storePlot(app.analyze.ensemble.endFitPlot);
            backup.analyze.ensemble.startLabel = storeText(app.analyze.ensemble.startLabel);
            backup.analyze.ensemble.endLabel = storeText(app.analyze.ensemble.endLabel);
            backup.analyze.ensemble.stepLabel = storeText(app.analyze.ensemble.stepLabel);
            
            backup.analyze.individual.pickChngptBtn = storeButton(app.analyze.individual.pickChngptBtn);
            backup.analyze.individual.resetBtn = storeButton(app.analyze.individual.resetBtn);
            backup.analyze.individual.snap = storeButton(app.analyze.individual.snap);
            backup.analyze.individual.snap.Value = app.analyze.individual.snap.Value;
            backup.analyze.individual.likelihoodAxes = storeAxes(app.analyze.individual.likelihoodAxes);
            backup.analyze.individual.surfPlot = storePlot(app.analyze.individual.surfPlot);
            
            backup.analyze.individual.defXonPlot = storePlot(app.analyze.individual.defXonPlot);
            backup.analyze.individual.defXoffPlot = storePlot(app.analyze.individual.defXoffPlot);
            backup.analyze.individual.defYonPlot = storePlot(app.analyze.individual.defYonPlot);
            backup.analyze.individual.defYoffPlot = storePlot(app.analyze.individual.defYoffPlot);
            backup.analyze.individual.XonPlot = storePlot(app.analyze.individual.XonPlot);
            backup.analyze.individual.XoffPlot = storePlot(app.analyze.individual.XoffPlot);
            backup.analyze.individual.YonPlot = storePlot(app.analyze.individual.YonPlot);
            backup.analyze.individual.YoffPlot = storePlot(app.analyze.individual.YoffPlot);
            backup.misc.simLines = storePlot(app.misc.simLines);
            backup.data.defaultEvents = app.data.defaultEvents;
            backup.data.windowStarts = app.data.windowStarts;
            backup.data.windowStops = app.data.windowStops;
            backup.data.defaultWindowStarts = app.data.defaultWindowStarts;
            backup.data.defaultWindowStops = app.data.defaultWindowStops;
            
            backup.analyze.individual.eventAxes = storeAxes(app.analyze.individual.eventAxes);
            backup.analyze.individual.eventAxes.YLabel.String = app.analyze.individual.eventAxes.YLabel.String;
            lh = getappdata(app.analyze.individual.eventAxes, 'listener');
            backup.analyze.individual.eventAxes.link = ~isempty(lh) && isvalid(lh);
            backup.analyze.individual.eventPlot = storePlot(app.analyze.individual.eventPlot);
            backup.analyze.individual.defTonPlot = storePlot(app.analyze.individual.defTonPlot);
            backup.analyze.individual.defToffPlot = storePlot(app.analyze.individual.defToffPlot);
            backup.analyze.individual.TonPlot = storePlot(app.analyze.individual.TonPlot);
            backup.analyze.individual.ToffPlot = storePlot(app.analyze.individual.ToffPlot);
            backup.analyze.individual.first = app.analyze.individual.first;
            backup.analyze.individual.T1 = app.analyze.individual.T1;
            backup.analyze.individual.defMean1Plot = storePlot(app.analyze.individual.defMean1Plot);
            backup.analyze.individual.defUpLim1Plot = storePlot(app.analyze.individual.defUpLim1Plot);
            backup.analyze.individual.defLowLim1Plot = storePlot(app.analyze.individual.defLowLim1Plot);
            backup.analyze.individual.defMean2Plot = storePlot(app.analyze.individual.defMean2Plot);
            backup.analyze.individual.defUpLim2Plot = storePlot(app.analyze.individual.defUpLim2Plot);
            backup.analyze.individual.defLowLim2Plot = storePlot(app.analyze.individual.defLowLim2Plot);
            backup.analyze.individual.defMean3Plot = storePlot(app.analyze.individual.defMean3Plot);
            backup.analyze.individual.defUpLim3Plot = storePlot(app.analyze.individual.defUpLim3Plot);
            backup.analyze.individual.defLowLim3Plot = storePlot(app.analyze.individual.defLowLim3Plot);
            backup.analyze.individual.mean1Plot = storePlot(app.analyze.individual.mean1Plot);
            backup.analyze.individual.upLim1Plot = storePlot(app.analyze.individual.upLim1Plot);
            backup.analyze.individual.lowLim1Plot = storePlot(app.analyze.individual.lowLim1Plot);
            backup.analyze.individual.mean2Plot = storePlot(app.analyze.individual.mean2Plot);
            backup.analyze.individual.upLim2Plot = storePlot(app.analyze.individual.upLim2Plot);
            backup.analyze.individual.lowLim2Plot = storePlot(app.analyze.individual.lowLim2Plot);
            backup.analyze.individual.mean3Plot = storePlot(app.analyze.individual.mean3Plot);
            backup.analyze.individual.upLim3Plot = storePlot(app.analyze.individual.upLim3Plot);
            backup.analyze.individual.lowLim3Plot = storePlot(app.analyze.individual.lowLim3Plot);
            
            backup.analyze.individual.panel.Visible = app.analyze.individual.panel.Visible;
            backup.analyze.individual.showRemovedBtn = storeButton(app.analyze.individual.showRemovedBtn);
            backup.analyze.individual.showRemovedBtn.Value = app.analyze.individual.showRemovedBtn.Value;
            backup.analyze.individual.currentEvent = app.analyze.individual.currentEvent;
            backup.analyze.individual.storeCurrentEvent = app.analyze.individual.storeCurrentEvent;
            backup.analyze.individual.currentEventTxt = storeText(app.analyze.individual.currentEventTxt);
            backup.analyze.individual.currentEventEdit = storeEdit(app.analyze.individual.currentEventEdit);
            backup.analyze.individual.ofNumTxt = storeText(app.analyze.individual.ofNumTxt);
            backup.analyze.individual.prevEventBtn = storeButton(app.analyze.individual.prevEventBtn);
            backup.analyze.individual.nextEventBtn = storeButton(app.analyze.individual.nextEventBtn);
            backup.analyze.individual.removeBtn = storeButton(app.analyze.individual.removeBtn);
            backup.analyze.individual.removeBtn.Value = app.analyze.individual.removeBtn.Value;
            backup.analyze.individual.label = storeText(app.analyze.individual.label);
            backup.analyze.individual.removedByDur = app.analyze.individual.removedByDur;
            backup.analyze.individual.removedBySep = app.analyze.individual.removedBySep;
            backup.analyze.individual.removedByUser = app.analyze.individual.removedByUser;
            backup.analyze.individual.includedByUser = app.analyze.individual.includedByUser;
            
            backup.controls.generateCovBtn = storeButton(app.controls.generateCovBtn);
            backup.controls.findEventsBtn = storeButton(app.controls.findEventsBtn);
            backup.controls.findEventsBtn.String = app.controls.findEventsBtn.String;
            backup.controls.useChangepointBtn = storeButton(app.controls.useChangepointBtn);
            backup.controls.useChangepointBtn.Value = app.controls.useChangepointBtn.Value;
            backup.controls.toggleEventsBtn = storeButton(app.controls.toggleEventsBtn);
            backup.controls.toggleEventsBtn.Value = app.controls.toggleEventsBtn.Value;
            backup.controls.flipDataBtn = storeButton(app.controls.flipDataBtn);
            backup.controls.flipDataBtn.Value = app.controls.flipDataBtn.Value;
            backup.controls.peakToPeakBtn = storeButton(app.controls.peakToPeakBtn);
            backup.controls.peakToPeakBtn.Value = app.controls.peakToPeakBtn.Value;
            backup.controls.thresholdBtn = storeButton(app.controls.thresholdBtn);
            backup.controls.automaticBtn = storeButton(app.controls.automaticBtn);
            backup.controls.automaticBtn.Value = app.controls.automaticBtn.Value;
            backup.controls.manualBtn = storeButton(app.controls.manualBtn);
            backup.controls.setManualBtn = storeButton(app.controls.setManualBtn);
            backup.controls.peak1Txt = storeText(app.controls.peak1Txt);
            backup.controls.peak1Edit = storeEdit(app.controls.peak1Edit);
            backup.controls.peak2Txt = storeText(app.controls.peak2Txt);
            backup.controls.peak2Edit = storeEdit(app.controls.peak2Edit);
            backup.controls.minTxt = storeText(app.controls.minTxt);
            backup.controls.minEdit = storeEdit(app.controls.minEdit);
            backup.controls.beadABtn = storeButton(app.controls.beadABtn);
            backup.controls.beadBBtn = storeButton(app.controls.beadBBtn);
            backup.controls.bothBeadsBtn = storeButton(app.controls.bothBeadsBtn);
            backup.controls.whichBeadBtnGrp.SelectedObject = app.controls.whichBeadBtnGrp.SelectedObject;
            
            backup.controls.deselect.panel.Visible = app.controls.deselect.panel.Visible;
            backup.controls.deselect.xpos = app.controls.deselect.xpos;
            backup.controls.deselect.linePlot = storePlot(app.controls.deselect.linePlot);
            backup.controls.deselect.leftPlot.Color = app.controls.deselect.leftPlot.Color;
            backup.controls.deselect.rightPlot.Color = app.controls.deselect.rightPlot.Color;
            backup.controls.deselect.patches = storePatch(app.controls.deselect.patches);
            backup.controls.deselect.first = app.controls.deselect.first;
            
            backup.controls.deselect.deselectBtn = storeButton(app.controls.deselect.deselectBtn);
            backup.controls.deselect.deselectBtn.Value = app.controls.deselect.deselectBtn.Value;
            backup.controls.deselect.deselectBtn.Callback = app.controls.deselect.deselectBtn.Callback;
            backup.controls.deselect.undoBtn = storeButton(app.controls.deselect.undoBtn);
            backup.controls.deselect.removeBtn = storeButton(app.controls.deselect.removeBtn);
            backup.controls.deselect.removeBtn.Callback = app.controls.deselect.removeBtn.Callback;
            backup.controls.deselect.saveBtn = storeButton(app.controls.deselect.saveBtn);
            backup.controls.deselect.indices = app.controls.deselect.indices;
            
            backup.controls.windows.panel.Visible = app.controls.windows.panel.Visible;
            backup.controls.windows.covwindow = app.controls.windows.covwindow;
            backup.controls.windows.covsmooth = app.controls.windows.covsmooth;
            backup.controls.windows.defcovwindow = app.controls.windows.defcovwindow;
            backup.controls.windows.defcovsmooth = app.controls.windows.defcovsmooth;
            backup.controls.windows.covwindowTxt = storeText(app.controls.windows.covwindowTxt);
            backup.controls.windows.covwindowEdit = storeText(app.controls.windows.covwindowEdit);
            backup.controls.windows.covsmoothTxt = storeText(app.controls.windows.covsmoothTxt);
            backup.controls.windows.covsmoothEdit = storeText(app.controls.windows.covsmoothEdit);
            backup.controls.windows.interactBtn = storeButton(app.controls.windows.interactBtn);
            backup.controls.windows.pt = storePlot(app.controls.windows.pt);
            backup.controls.windows.resetBtn = storeButton(app.controls.windows.resetBtn);
            
            backup.controls.cutoffs.defMinSep = app.controls.cutoffs.defMinSep;
            backup.controls.cutoffs.defMinDur = app.controls.cutoffs.defMinDur;
            backup.controls.cutoffs.minDurTxt = storeText(app.controls.cutoffs.minDurTxt);
            backup.controls.cutoffs.minDurEdit = storeEdit(app.controls.cutoffs.minDurEdit);
            backup.controls.cutoffs.minSepTxt = storeText(app.controls.cutoffs.minSepTxt);
            backup.controls.cutoffs.minSepEdit = storeEdit(app.controls.cutoffs.minSepEdit);
            backup.controls.cutoffs.resetBtn = storeButton(app.controls.cutoffs.resetBtn);
            
            backup.misc.active.allEvents.cov = app.misc.active.allEvents.cov;
            backup.misc.active.allEvents.peak1 = app.misc.active.allEvents.peak1;
            backup.misc.active.allEvents.peak2 = app.misc.active.allEvents.peak2;
            backup.misc.active.allEvents.min = app.misc.active.allEvents.min;
            backup.misc.active.allEvents.whichBeads = app.misc.active.allEvents.whichBeads;
            backup.misc.active.allEvents.useChangepoint = app.misc.active.allEvents.useChangepoint;
            
            backup.excel.listbox.String = app.excel.listbox.String;
            backup.excel.listbox.Value = app.excel.listbox.Value;
            backup.excel.listbox.Enable = app.excel.listbox.Enable;
            backup.excel.files = app.excel.files;
            backup.excel.createBtn = storeButton(app.excel.createBtn);
            backup.excel.findBtn = storeButton(app.excel.findBtn);
            backup.excel.addBtn = storeButton(app.excel.addBtn);
            
            backup.misc.toBeEnabled = app.misc.toBeEnabled;
            backup.misc.disabled = app.misc.disabled;
            backup.misc.skip = app.misc.skip;
            enableable = app.getEnableable;
            for e = 1:numel(enableable)
                backup.enableable(e).Enable = enableable(e).Enable;
            end
            backup.load.currentTask = storeText(app.load.currentTask);
            backup.load2.currentTask = storeText(app.load2.currentTask);
            if ~app.misc.axToolbars
                backup.misc.menus.Enable = app.misc.menus(1).Enable;
            end
            
            backup.allEvents = app.allEvents;
            backup.activeEvents = app.activeEvents;
            backup.minSep = app.minSep;
            backup.minDur = app.minDur;
            backup.autop1 = app.autop1;
            backup.autop2 = app.autop2;
            backup.automin = app.automin;
            backup.manp1 = app.manp1;
            backup.manp2 = app.manp2;
            backup.manmin = app.manmin;
            backup.useman = app.useman;
            backup.usemin = app.usemin;
            backup.flipped = app.flipped;
            backup.showEvents = app.showEvents;
            backup.showRemoved = app.showRemoved;
            
            function copyB = storeButton(b)
                copyB.Visible = b.Visible;
                copyB.Enable = b.Enable;
            end
            
            function copyT = storeText(t)
                copyT.Visible = t.Visible;
                if isprop(t, 'Enable'), copyT.Enable = t.Enable; end
                copyT.String = t.String;
                copyT.Position = t.Position;
            end
            
            function copyE = storeEdit(e)
                copyE.Visible = e.Visible;
                copyE.Enable = e.Enable;
                copyE.String = e.String;
            end
            
            function copyA = storeAxes(a)
                copyA.Visible = a.Visible;
                copyA.XColor = a.XColor;
                copyA.YColor = a.YColor;
                copyA.ZColor = a.ZColor;
                copyA.Position = a.Position;
                copyA.View = a.View;
                copyA.Limits = axis(a);
                copyA.XLimMode = a.XLimMode;
                copyA.YLimMode = a.YLimMode;
                copyA.ZLimMode = a.ZLimMode;
                copyA.appdata = getappdata(a);
                if app.misc.axToolbars
                    copyA.Toolbar.Visible = a.Toolbar.Visible;
                end
            end
            
            function copyP = storePlot(p)
                if isempty(p), copyP = []; return, end
                copyP(numel(p)) = struct;
                for i = 1:numel(p)
                    copyP(i).Visible = p(i).Visible;
                    copyP(i).XData = p(i).XData;
                    copyP(i).YData = p(i).YData;
                    if isprop(p(i), 'ZData'), copyP(i).ZData = p(i).ZData; end
                    copyP(i).Parent = p(i).Parent;
                end
            end
            
            function copyP = storePatch(p)
                if isempty(p), copyP = []; return, end
                copyP(numel(p)) = struct;
                for i = 1:numel(p)
                    copyP(i).Visible = p(i).Visible;
                    copyP(i).XData = p(i).XData;
                    copyP(i).YData = p(i).YData;
                    copyP(i).Parent = p(i).Parent;
                end
            end
        end
        
        function restoreBackupTab1(app, backup, mode)
            % Restores tab 1 properties to the values contained in the struct backup. The
            % subset of properties which are updated depends on mode.
            
            app.misc.postListener.Enabled = false;
            
            if nargin < 3
                mode = 'all';
            end
            
            figQ = any(strcmp(mode, {'all', 'fullReset'}));
            loadQ = any(strcmp(mode, {'all', 'fullReset'}));
            mainQ = any(strcmp(mode, {'all', 'fullReset'}));
            covQ = any(strcmp(mode, {'all', 'fullReset'}));
            covHistQ = any(strcmp(mode, {'all', 'fullReset'}));
            miscQ = any(strcmp(mode, {'all', 'fullReset', 'newEventsReset', 'noActiveEventsReset'}));
            ensembleQ = any(strcmp(mode, {'all', 'fullReset', 'newEventsReset', 'noActiveEventsReset'}));
            individualQ = any(strcmp(mode, {'all', 'fullReset', 'newEventsReset'}));
            controlsQ = any(strcmp(mode, {'all', 'fullReset'}));
            deselectQ = any(strcmp(mode, {'all', 'fullReset'}));
            eventsQ = any(strcmp(mode, {'all', 'fullReset', 'newEventsReset'}));
            excelQ = any(strcmp(mode, {'all', 'fullReset'}));
            excel_filesQ = any(strcmp(mode, 'all'));
            dataQ = any(strcmp(mode, {'all', 'fullReset'}));
            generalQ = any(strcmp(mode, {'all'}));
            
            if figQ
                app.misc.fig.WindowButtonMotionFcn = backup.misc.fig.WindowButtonMotionFcn;
                app.misc.fig.WindowButtonDownFcn = backup.misc.fig.WindowButtonDownFcn;
                app.misc.fig.WindowKeyPressFcn = backup.misc.fig.WindowKeyPressFcn;
                app.misc.fig.SizeChangedFcn = backup.misc.fig.SizeChangedFcn;
                app.misc.storeClickCallback = backup.misc.storeClickCallback;
                app.misc.storeMotionCallback = backup.misc.storeMotionCallback;
                app.misc.fig.Pointer = backup.misc.fig.Pointer;
            end
            
            if loadQ
                restoreButton(app.load.inputBtn, backup.load.inputBtn);
                app.load.path = backup.load.path;
                app.load.name = backup.load.name;
                app.load.ext = backup.load.ext;
                app.load.header = backup.load.header;
                app.data.origData = backup.data.origData;
                app.data.origBeads = backup.data.origBeads;
                app.data.time = backup.data.time;
                app.load.CALA = backup.load.CALA;
                restoreText(app.load.CALATxt, backup.load.CALATxt);
                app.load.CALB = backup.load.CALB;
                restoreText(app.load.CALBTxt, backup.load.CALBTxt);
                app.load.KA = backup.load.KA;
                restoreText(app.load.KATxt, backup.load.KATxt);
                app.load.KB = backup.load.KB;
                restoreText(app.load.KBTxt, backup.load.KBTxt);
                app.data.fs = backup.data.fs;
                restoreText(app.load.fsTxt, backup.load.fsTxt);
                restoreText(app.load.filenameTxt, backup.load.filenameTxt);
                restoreButton(app.load.skipBtn, backup.load.skipBtn);
                app.load.skipBtn.Position = backup.load.skipBtn.Position;
                app.load.taskMaxWidth = backup.load.taskMaxWidth;
                app.misc.simInput = backup.misc.simInput;
                app.misc.realEvents = backup.misc.realEvents;
            end
            
            if mainQ
                restoreAxes(app.load.main.mainAxes, backup.load.main.mainAxes);
                app.load.main.mainAxes.XLabel.Visible = backup.load.main.mainAxes.XLabel.Visible;
                app.load.main.mainAxes.XTickMode = backup.load.main.mainAxes.XTickMode;
                app.load.main.mainAxes.XTickLabelMode = backup.load.main.mainAxes.XTickLabelMode;
                app.load.main.ylim = backup.load.main.ylim;
                restorePlot(app.load.main.bead1Plot, backup.load.main.bead1Plot);
                restorePlot(app.load.main.bead2Plot, backup.load.main.bead2Plot);
                if eventsQ
                    app.load.main.patches = restorePatch(app.load.main.patches, backup.load.main.patches,...
                        'FaceColor', app.misc.colors.mainPatches,...
                        'FaceAlpha', 0.3);
                    app.misc.simPatches = restorePatch(app.misc.simPatches, backup.misc.simPatches,...
                        'FaceColor', app.misc.colors.simPatches,...
                        'EdgeColor', 'none',...
                        'FaceAlpha', 0.3);
                end
                restoreText(app.load.main.label, backup.load.main.label);
                app.load.main.prevIn = backup.load.main.prevIn;
            end
            
            if covQ
                restoreAxes(app.load.cov.covAxes, backup.load.cov.covAxes);
                restorePlot(app.load.cov.covPlot, backup.load.cov.covPlot);
                restorePlot(app.load.cov.peak1Plot, backup.load.cov.peak1Plot);
                restorePlot(app.load.cov.peak2Plot, backup.load.cov.peak2Plot);
                restorePlot(app.load.cov.minPlot, backup.load.cov.minPlot);
                if eventsQ
                    app.load.cov.patches = restorePatch(app.load.cov.patches, backup.load.cov.patches,...
                        'FaceColor', app.misc.colors.covPatches,...
                        'FaceAlpha', 0.3);
                end
                app.load.cov.cov = backup.load.cov.cov;
            end
            
            if covHistQ
                restoreAxes(app.analyze.covHist.covHistAxes, backup.analyze.covHist.covHistAxes);
                app.analyze.covHist.covHistPlot.Data = backup.analyze.covHist.covHistPlot.Data;
                app.analyze.covHist.covHistPlot.Visible = backup.analyze.covHist.covHistPlot.Visible;
                app.analyze.covHist.covHistPlot.NumBins = backup.analyze.covHist.covHistPlot.NumBins;
                restorePlot(app.analyze.covHist.peak1Plot, backup.analyze.covHist.peak1Plot);
                restorePlot(app.analyze.covHist.peak2Plot, backup.analyze.covHist.peak2Plot);
                restorePlot(app.analyze.covHist.minPlot, backup.analyze.covHist.minPlot);
            end
            
            if miscQ
                restoreAxes(app.analyze.misc.miscAxes, backup.analyze.misc.miscAxes);
                app.analyze.misc.miscAxes.YScale = backup.analyze.misc.miscAxes.YScale;
                app.analyze.misc.miscAxes.YLabel.String = backup.analyze.misc.miscAxes.YLabel.String;
                app.analyze.misc.miscAxes.XLabel.String = backup.analyze.misc.miscAxes.XLabel.String;
                restorePlot(app.analyze.misc.miscPlot, backup.analyze.misc.miscPlot);
                app.analyze.misc.miscPlot.Marker = backup.analyze.misc.miscPlot.Marker;
                app.analyze.misc.miscPlot.LineStyle = backup.analyze.misc.miscPlot.LineStyle;
                restorePlot(app.analyze.misc.realPlot, backup.analyze.misc.realPlot);
                restorePlot(app.analyze.misc.miscFitPlot, backup.analyze.misc.miscFitPlot);
                restoreText(app.analyze.misc.miscLabel, backup.analyze.misc.miscLabel);
                app.analyze.misc.data = backup.analyze.misc.data;
                app.analyze.misc.kDur = backup.analyze.misc.kDur;
                app.analyze.misc.meanStep1 = backup.analyze.misc.meanStep1;
                app.analyze.misc.stdStep1 = backup.analyze.misc.stdStep1;
                app.analyze.misc.meanStep2 = backup.analyze.misc.meanStep2;
                app.analyze.misc.stdStep2 = backup.analyze.misc.stdStep2;
                app.analyze.misc.meanTotalStep = backup.analyze.misc.meanTotalStep;
                app.analyze.misc.stdTotalStep = backup.analyze.misc.stdTotalStep;
            end
            
            if ensembleQ
                restoreAxes(app.analyze.ensemble.ensembleAxes, backup.analyze.ensemble.ensembleAxes);
                app.analyze.ensemble.ensembleAxes.YLabel.String = backup.analyze.ensemble.ensembleAxes.YLabel.String;
                app.analyze.ensemble.globalFit = backup.analyze.ensemble.globalFit;
                app.analyze.ensemble.start = backup.analyze.ensemble.start;
                app.analyze.ensemble.end = backup.analyze.ensemble.end;
                app.analyze.ensemble.startFit = backup.analyze.ensemble.startFit;
                app.analyze.ensemble.endFit = backup.analyze.ensemble.endFit;
                app.analyze.ensemble.p = backup.analyze.ensemble.p;
                restorePlot(app.analyze.ensemble.startPlot, backup.analyze.ensemble.startPlot);
                restorePlot(app.analyze.ensemble.endPlot, backup.analyze.ensemble.endPlot);
                restorePlot(app.analyze.ensemble.startFitPlot, backup.analyze.ensemble.startFitPlot);
                restorePlot(app.analyze.ensemble.endFitPlot, backup.analyze.ensemble.endFitPlot);
                restoreText(app.analyze.ensemble.startLabel, backup.analyze.ensemble.startLabel);
                restoreText(app.analyze.ensemble.endLabel, backup.analyze.ensemble.endLabel);
                restoreText(app.analyze.ensemble.stepLabel, backup.analyze.ensemble.stepLabel);
            end
            
            if individualQ
                restoreButton(app.analyze.individual.pickChngptBtn, backup.analyze.individual.pickChngptBtn);
                restoreButton(app.analyze.individual.resetBtn, backup.analyze.individual.resetBtn);
                restoreButton(app.analyze.individual.snap, backup.analyze.individual.snap);
                app.analyze.individual.snap.Value = backup.analyze.individual.snap.Value;
                restoreAxes(app.analyze.individual.likelihoodAxes, backup.analyze.individual.likelihoodAxes);
                restorePlot(app.analyze.individual.surfPlot, backup.analyze.individual.surfPlot);
                restorePlot(app.analyze.individual.defXonPlot, backup.analyze.individual.defXonPlot);
                restorePlot(app.analyze.individual.defXoffPlot, backup.analyze.individual.defXoffPlot);
                restorePlot(app.analyze.individual.defYonPlot, backup.analyze.individual.defYonPlot);
                restorePlot(app.analyze.individual.defYoffPlot, backup.analyze.individual.defYoffPlot);
                restorePlot(app.analyze.individual.XonPlot, backup.analyze.individual.XonPlot);
                restorePlot(app.analyze.individual.XoffPlot, backup.analyze.individual.XoffPlot);
                restorePlot(app.analyze.individual.YonPlot, backup.analyze.individual.YonPlot);
                restorePlot(app.analyze.individual.YoffPlot, backup.analyze.individual.YoffPlot);
                app.misc.simLines = restorePlot(app.misc.simLines, backup.misc.simLines,...
                    'Color', app.misc.colors.simLines);
                app.data.defaultEvents = backup.data.defaultEvents;
                app.data.windowStarts = backup.data.windowStarts;
                app.data.windowStops = backup.data.windowStops;
                app.data.defaultWindowStarts = backup.data.defaultWindowStarts;
                app.data.defaultWindowStops = backup.data.defaultWindowStops;
                
                delete(getappdata(app.analyze.individual.eventAxes, 'listener'));
                restoreAxes(app.analyze.individual.eventAxes, backup.analyze.individual.eventAxes);
                app.analyze.individual.eventAxes.YLabel.String = backup.analyze.individual.eventAxes.YLabel.String;
                restorePlot(app.analyze.individual.eventPlot, backup.analyze.individual.eventPlot);
                restorePlot(app.analyze.individual.defTonPlot, backup.analyze.individual.defTonPlot);
                restorePlot(app.analyze.individual.defToffPlot, backup.analyze.individual.defToffPlot);
                restorePlot(app.analyze.individual.TonPlot, backup.analyze.individual.TonPlot);
                restorePlot(app.analyze.individual.ToffPlot, backup.analyze.individual.ToffPlot);
                app.analyze.individual.first = backup.analyze.individual.first;
                app.analyze.individual.T1 = backup.analyze.individual.T1;
                restorePlot(app.analyze.individual.defMean1Plot, backup.analyze.individual.defMean1Plot);
                restorePlot(app.analyze.individual.defUpLim1Plot, backup.analyze.individual.defUpLim1Plot);
                restorePlot(app.analyze.individual.defLowLim1Plot, backup.analyze.individual.defLowLim1Plot);
                restorePlot(app.analyze.individual.defMean2Plot, backup.analyze.individual.defMean2Plot);
                restorePlot(app.analyze.individual.defUpLim2Plot, backup.analyze.individual.defUpLim2Plot);
                restorePlot(app.analyze.individual.defLowLim2Plot, backup.analyze.individual.defLowLim2Plot);
                restorePlot(app.analyze.individual.defMean3Plot, backup.analyze.individual.defMean3Plot);
                restorePlot(app.analyze.individual.defUpLim3Plot, backup.analyze.individual.defUpLim3Plot);
                restorePlot(app.analyze.individual.defLowLim3Plot, backup.analyze.individual.defLowLim3Plot);
                restorePlot(app.analyze.individual.mean1Plot, backup.analyze.individual.mean1Plot);
                restorePlot(app.analyze.individual.upLim1Plot, backup.analyze.individual.upLim1Plot);
                restorePlot(app.analyze.individual.lowLim1Plot, backup.analyze.individual.lowLim1Plot);
                restorePlot(app.analyze.individual.mean2Plot, backup.analyze.individual.mean2Plot);
                restorePlot(app.analyze.individual.upLim2Plot, backup.analyze.individual.upLim2Plot);
                restorePlot(app.analyze.individual.lowLim2Plot, backup.analyze.individual.lowLim2Plot);
                restorePlot(app.analyze.individual.mean3Plot, backup.analyze.individual.mean3Plot);
                restorePlot(app.analyze.individual.upLim3Plot, backup.analyze.individual.upLim3Plot);
                restorePlot(app.analyze.individual.lowLim3Plot, backup.analyze.individual.lowLim3Plot);
                
                app.analyze.individual.panel.Visible = backup.analyze.individual.panel.Visible;
                restoreButton(app.analyze.individual.showRemovedBtn, backup.analyze.individual.showRemovedBtn);
                app.analyze.individual.showRemovedBtn.Value = backup.analyze.individual.showRemovedBtn.Value;
                app.analyze.individual.currentEvent = backup.analyze.individual.currentEvent;
                app.analyze.individual.storeCurrentEvent = backup.analyze.individual.storeCurrentEvent;
                restoreText(app.analyze.individual.currentEventTxt, backup.analyze.individual.currentEventTxt);
                restoreEdit(app.analyze.individual.currentEventEdit, backup.analyze.individual.currentEventEdit);
                restoreText(app.analyze.individual.ofNumTxt, backup.analyze.individual.ofNumTxt);
                restoreButton(app.analyze.individual.prevEventBtn, backup.analyze.individual.prevEventBtn);
                restoreButton(app.analyze.individual.nextEventBtn, backup.analyze.individual.nextEventBtn);
                restoreButton(app.analyze.individual.removeBtn, backup.analyze.individual.removeBtn);
                app.analyze.individual.removeBtn.Value = backup.analyze.individual.removeBtn.Value;
                restoreText(app.analyze.individual.label, backup.analyze.individual.label);
                app.analyze.individual.removedByDur = backup.analyze.individual.removedByDur;
                app.analyze.individual.removedBySep = backup.analyze.individual.removedBySep;
                app.analyze.individual.removedByUser = backup.analyze.individual.removedByUser;
                app.analyze.individual.includedByUser = backup.analyze.individual.includedByUser;
            end
            
            if controlsQ
                restoreButton(app.controls.generateCovBtn, backup.controls.generateCovBtn);
                restoreButton(app.controls.findEventsBtn, backup.controls.findEventsBtn);
                app.controls.findEventsBtn.String = backup.controls.findEventsBtn.String;
                restoreButton(app.controls.useChangepointBtn, backup.controls.useChangepointBtn);
                app.controls.useChangepointBtn.Value = backup.controls.useChangepointBtn.Value;
                restoreButton(app.controls.toggleEventsBtn, backup.controls.toggleEventsBtn);
                app.controls.toggleEventsBtn.Value = backup.controls.toggleEventsBtn.Value;
                restoreButton(app.controls.flipDataBtn, backup.controls.flipDataBtn);
                app.controls.flipDataBtn.Value = backup.controls.flipDataBtn.Value;
                restoreButton(app.controls.peakToPeakBtn, backup.controls.peakToPeakBtn);
                app.controls.peakToPeakBtn.Value = backup.controls.peakToPeakBtn.Value;
                restoreButton(app.controls.thresholdBtn, backup.controls.thresholdBtn);
                restoreButton(app.controls.automaticBtn, backup.controls.automaticBtn);
                app.controls.automaticBtn.Value = backup.controls.automaticBtn.Value;
                restoreButton(app.controls.manualBtn, backup.controls.manualBtn);
                restoreButton(app.controls.setManualBtn, backup.controls.setManualBtn);
                
                restoreText(app.controls.peak1Txt, backup.controls.peak1Txt);
                restoreEdit(app.controls.peak1Edit, backup.controls.peak1Edit);
                restoreText(app.controls.peak2Txt, backup.controls.peak2Txt);
                restoreEdit(app.controls.peak2Edit, backup.controls.peak2Edit);
                restoreText(app.controls.minTxt, backup.controls.minTxt);
                restoreEdit(app.controls.minEdit, backup.controls.minEdit);
                
                restoreButton(app.controls.beadABtn, backup.controls.beadABtn);
                restoreButton(app.controls.beadBBtn, backup.controls.beadBBtn);
                restoreButton(app.controls.bothBeadsBtn, backup.controls.bothBeadsBtn);
                app.controls.whichBeadBtnGrp.SelectedObject = backup.controls.whichBeadBtnGrp.SelectedObject;
                
                app.controls.windows.panel.Visible = backup.controls.windows.panel.Visible;
                app.controls.windows.covwindow = backup.controls.windows.covwindow;
                app.controls.windows.covsmooth = backup.controls.windows.covsmooth;
                app.controls.windows.defcovwindow = backup.controls.windows.defcovwindow;
                app.controls.windows.defcovsmooth = backup.controls.windows.defcovsmooth;
                restoreText(app.controls.windows.covwindowTxt, backup.controls.windows.covwindowTxt);
                restoreText(app.controls.windows.covwindowEdit, backup.controls.windows.covwindowEdit);
                restoreText(app.controls.windows.covsmoothTxt, backup.controls.windows.covsmoothTxt);
                restoreText(app.controls.windows.covsmoothEdit, backup.controls.windows.covsmoothEdit);
                restoreButton(app.controls.windows.interactBtn, backup.controls.windows.interactBtn);
                restorePlot(app.controls.windows.pt, backup.controls.windows.pt);
                restoreButton(app.controls.windows.resetBtn, backup.controls.windows.resetBtn);
                
                app.controls.cutoffs.defMinSep = backup.controls.cutoffs.defMinSep;
                app.controls.cutoffs.defMinDur = backup.controls.cutoffs.defMinDur;
                restoreText(app.controls.cutoffs.minDurTxt, backup.controls.cutoffs.minDurTxt);
                restoreEdit(app.controls.cutoffs.minDurEdit, backup.controls.cutoffs.minDurEdit);
                restoreText(app.controls.cutoffs.minSepTxt, backup.controls.cutoffs.minSepTxt);
                restoreEdit(app.controls.cutoffs.minSepEdit, backup.controls.cutoffs.minSepEdit);
                restoreButton(app.controls.cutoffs.resetBtn, backup.controls.cutoffs.resetBtn);
            end
            
            if deselectQ
                app.controls.deselect.panel.Visible = backup.controls.deselect.panel.Visible;
                app.controls.deselect.xpos = backup.controls.deselect.xpos;
                restorePlot(app.controls.deselect.linePlot, backup.controls.deselect.linePlot);
                app.controls.deselect.leftPlot.Color = backup.controls.deselect.leftPlot.Color;
                app.controls.deselect.rightPlot.Color = backup.controls.deselect.rightPlot.Color;
                app.controls.deselect.patches = restorePatch(app.controls.deselect.patches, backup.controls.deselect.patches,...
                    'FaceColor', app.misc.colors.deselect,...
                    'FaceAlpha', 0.3);
                app.controls.deselect.first = backup.controls.deselect.first;
                
                restoreButton(app.controls.deselect.deselectBtn, backup.controls.deselect.deselectBtn);
                app.controls.deselect.deselectBtn.Value = backup.controls.deselect.deselectBtn.Value;
                app.controls.deselect.deselectBtn.Callback = backup.controls.deselect.deselectBtn.Callback;
                restoreButton(app.controls.deselect.undoBtn, backup.controls.deselect.undoBtn);
                restoreButton(app.controls.deselect.removeBtn, backup.controls.deselect.removeBtn);
                app.controls.deselect.removeBtn.Callback = backup.controls.deselect.removeBtn.Callback;
                restoreButton(app.controls.deselect.saveBtn, backup.controls.deselect.saveBtn);
                app.controls.deselect.indices = backup.controls.deselect.indices;
            end
            
            if eventsQ
                app.misc.active.allEvents.cov = backup.misc.active.allEvents.cov;
                app.misc.active.allEvents.peak1 = backup.misc.active.allEvents.peak1;
                app.misc.active.allEvents.peak2 = backup.misc.active.allEvents.peak2;
                app.misc.active.allEvents.min = backup.misc.active.allEvents.min;
                app.misc.active.allEvents.whichBeads = backup.misc.active.allEvents.whichBeads;
                app.misc.active.allEvents.useChangepoint = backup.misc.active.allEvents.useChangepoint;
            end
            
            if mainQ
                restoreAxesLimits(app.load.main.mainAxes, backup.load.main.mainAxes);
            end
            if covQ
                restoreAxesLimits(app.load.cov.covAxes, backup.load.cov.covAxes);
            end
            if covHistQ
                restoreAxesLimits(app.analyze.covHist.covHistAxes, backup.analyze.covHist.covHistAxes);
            end
            if miscQ
                restoreAxesLimits(app.analyze.misc.miscAxes, backup.analyze.misc.miscAxes);
            end
            if ensembleQ
                restoreAxesLimits(app.analyze.ensemble.ensembleAxes, backup.analyze.ensemble.ensembleAxes);
            end
            if individualQ
                restoreAxesLimits(app.analyze.individual.likelihoodAxes, backup.analyze.individual.likelihoodAxes);
                restoreAxesLimits(app.analyze.individual.eventAxes, backup.analyze.individual.eventAxes);
                if backup.analyze.individual.eventAxes.link
                    lh = addlistener(app.analyze.individual.eventAxes, 'XLim', 'PostSet', @app.eventXAxisUpdate);
                    setappdata(app.analyze.individual.eventAxes, 'listener', lh);
                end
            end
            
            if excelQ
                app.excel.listbox.Enable = backup.excel.listbox.Enable;
                restoreButton(app.excel.createBtn, backup.excel.createBtn);
                restoreButton(app.excel.findBtn, backup.excel.findBtn);
                restoreButton(app.excel.addBtn, backup.excel.addBtn);
                if excel_filesQ
                    app.excel.listbox.String = backup.excel.listbox.String;
                    app.excel.listbox.Value = backup.excel.listbox.Value;
                    app.excel.files = backup.excel.files;
                end
            end
            
            if eventsQ
                app.allEvents = backup.allEvents;
                app.activeEvents = backup.activeEvents;
                app.showEvents = backup.showEvents;
                app.showRemoved = backup.showRemoved;
            end
            
            if generalQ
                app.misc.fig.Position = backup.misc.fig.Position;
                app.misc.toBeEnabled = backup.misc.toBeEnabled;
                app.misc.disabled = backup.misc.disabled;
                app.misc.skip = backup.misc.skip;
                enableable = app.getEnableable;
                for e = 1:numel(enableable)
                    enableable(e).Enable = backup.enableable(e).Enable;
                end
                restoreText(app.load.currentTask, backup.load.currentTask);
                restoreText(app.load2.currentTask, backup.load2.currentTask);
                if ~app.misc.axToolbars
                    set(app.misc.menus, 'Enable', backup.misc.menus.Enable);
                end
            end
            
            if dataQ                
                app.minSep = backup.minSep;
                app.minDur = backup.minDur;
                app.autop1 = backup.autop1;
                app.autop2 = backup.autop2;
                app.automin = backup.automin;
                app.manp1 = backup.manp1;
                app.manp2 = backup.manp2;
                app.manmin = backup.manmin;
                app.useman = backup.useman;
                app.usemin = backup.usemin;
                app.flipped = backup.flipped;
            end
            
            app.misc.postListener.Enabled = true;
            drawnow
            
            function restoreButton(b, copyB)
                b.Visible = copyB.Visible;
                b.Enable = copyB.Enable;
            end
            
            function restoreText(t, copyT)
                t.Visible = copyT.Visible;
                if isprop(t, 'Enable'), t.Enable = copyT.Enable; end
                t.String = copyT.String;
                t.Position = copyT.Position;
            end
            
            function restoreEdit(e, copyE)
                e.Visible = copyE.Visible;
                e.Enable = copyE.Enable;
                e.String = copyE.String;
            end
            
            function restoreAxes(a, copyA)
                a.Visible = copyA.Visible;
                a.XColor = copyA.XColor;
                a.YColor = copyA.YColor;
                a.ZColor = copyA.ZColor;
                a.Position = copyA.Position;
                a.View = copyA.View;
                a.XLimMode = 'auto';
                a.YLimMode = 'auto';
                a.ZLimMode = 'auto';
                if app.misc.axToolbars
                    a.Toolbar.Visible = copyA.Toolbar.Visible;
                end
                % Restore appdata, which stores links, default limits, maybe other stuff:
                names = fieldnames(getappdata(a));
                for iName = 1:numel(names)
                    rmappdata(a, names{iName});
                end
                oldNames = fieldnames(copyA.appdata);
                for iName = 1:numel(oldNames)
                    setappdata(a, oldNames{iName}, copyA.appdata.(oldNames{iName}));
                end
            end
            
            function restoreAxesLimits(a, copyA)
                if strcmp(a.Visible, 'on')
                    axis(a, copyA.Limits);
                end
                a.XLimMode = copyA.XLimMode;
                a.YLimMode = copyA.YLimMode;
                a.ZLimMode = copyA.ZLimMode;
            end
            
            function p = restorePlot(p, copyP, varargin)
                if isempty(copyP), delete(p); p = gobjects(0); return, end
                if numel(p) < numel(copyP), delete(p); p = gobjects(numel(copyP),1); end
                if numel(p) > numel(copyP), delete(p(numel(copyP)+1:end)); p = p(1:numel(copyP)); end
                for i = 1:numel(copyP)
                    if ~ishghandle(p(i))
                        p(i) = plot(copyP(i).Parent, copyP(i).XData, copyP(i).YData,...
                            'Visible', copyP(i).Visible,...        
                            'HitTest', 'off',...
                            varargin{:});
                    else
                        p(i).Visible = copyP(i).Visible;
                        p(i).XData = copyP(i).XData;
                        p(i).YData = copyP(i).YData;
                        if isfield(copyP(i), 'ZData'), p(i).ZData = copyP(i).ZData; end
                    end
                end
            end
            
            function p = restorePatch(p, copyP, varargin)
                if isempty(copyP), delete(p); p = gobjects(0); return, end
                if numel(p) < numel(copyP), delete(p); p = gobjects(numel(copyP),1); end
                if numel(p) > numel(copyP), delete(p(numel(copyP)+1:end)); p = p(1:numel(copyP)); end
                for i = 1:numel(copyP)
                    if ~ishghandle(p(i))
                        p(i) = patch(copyP(i).Parent,...
                            'XData', copyP(i).XData,...
                            'YData', copyP(i).YData,...
                            'Visible', copyP(i).Visible,...
                            'HitTest', 'off',...
                            varargin{:});
                    else
                        p(i).Visible = copyP(i).Visible;
                        p(i).XData = copyP(i).XData;
                        p(i).YData = copyP(i).YData;
                    end
                end
            end
        end
        
        function fullReset(app)
            % Reset tab 1 properties to exactly as they were after initialize().
            
            app.restoreBackupTab1(app.misc.initBackupTab1, 'fullReset');
        end
        
        function newEventsReset(app)
            % Reset the tab 1 properties which get set after events are found.
            
            app.restoreBackupTab1(app.misc.initBackupTab1, 'newEventsReset');
        end
        
        function noActiveEventsReset(app)
            % Hide axes and plots in tab 1 which require there to be at least one active
            % event.
            
            app.restoreBackupTab1(app.misc.initBackupTab1, 'noActiveEventsReset');
        end
        
        function chooseTxtCallback(app, ~, ~)
            % Called when user clicks load.inputBtn. Allows user to pick a .txt file and
            % then updates load panel of tab 1.
            
            backup = app.storeBackupTab1;
            try
                % Let user pick .txt file.
                [name, path] = uigetfile('.txt');
                if name == 0
                    return
                end
                app.fullReset;
                store = app.disable;
                drawnow;
                
                % Read data.
                app.updateCurrentTask('Reading data.');
                [app.load.path, app.load.name, app.load.ext] = fileparts([path name]);
                if ~strcmp(app.load.ext, '.txt')
                    uiwait(errordlg('Input must be a .txt file.', 'Error Dialog', 'modal'))
                    app.restoreBackupTab1(backup);
                    return
                end
                warning('off', 'MATLAB:table:ModifiedAndSavedVarnames') % Suppress warning output.
                try
                    % Try reading data.
                    try
                        % Try using the data input function.
                        [app.data.fs, app.load.KA, app.load.KB,...
                            app.load.CALA, app.load.CALB,...
                            app.data.origData, app.load.header] = app.load.txtInputFcn([path name]);
                    catch ME
                        % Display any error to user and exit.
                        app.unexpectedError(ME, sprintf(['There was an error while reading your data:\n\n'...
                            '%s\n\nDouble check that your .txt file is formatted correctly.'], ME.message));
                        warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
                        app.restoreBackupTab1(backup);
                        return
                    end
                    % If function ran successfully, now test outputs.
                    testNumeric = @(a, s) assert(isa(a, 'numeric') && ~isnan(a),...
                        'chooseTxtCallback:outputFailed',...
                        'The %s returned by the data input function was not a number.',...
                        s);
                    testNumeric(app.load.KA, 'value of K for bead A');
                    testNumeric(app.load.KB, 'value of K for bead B');
                    testNumeric(app.load.CALA, 'value of CAL for bead A');
                    testNumeric(app.load.CALB, 'value of CAL for bead B');
                    testNumeric(app.data.fs, 'sample rate');
                    assert(istable(app.data.origData),...
                        'chooseTxtCallback:outputFailed',...
                        'The final variable returned by the data input function had type %s but should be a table.',...
                        class(app.data.origData));
                    i = ismember({'BeadAPos','BeadBPos'}, app.data.origData.Properties.VariableNames);
                    assert(i(1),...
                        'chooseTxtCallback:outputFailed',...
                        'The table returned by the data input function did not contain a column labeled "BeadAPos".');
                    assert(i(2),...
                        'chooseTxtCallback:outputFailed',...
                        'The table returned by the data input function did not contain a column labeled "BeadBPos".');
                    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
                catch ME
                    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
                    if strcmp(ME.identifier, 'chooseTxtCallback:outputFailed')
                        % If any of the outputs are invalid:
                        app.unexpectedError(ME, sprintf('%s Double check that your .txt file is formatted correctly.', ME.message));
                    else
                        % If there was some other error:
                        app.unexpectedError(ME, sprintf(['There was an unexpected error when '...
                            'trying to read your .txt file:\n\n%s'], ME.message));
                    end
                    app.restoreBackupTab1(backup);
                    return
                end

                app.misc.simInput = ismember('Key', app.data.origData.Properties.VariableNames);

                % Here, we know that the input was valid. Store data, display header info, and
                % plot beads.
                app.data.origBeads = [app.data.origData.BeadAPos*app.load.CALA,...
                    app.data.origData.BeadBPos*app.load.CALB];
                app.controls.deselect.indices = true(size(app.data.origBeads, 1), 1);
                app.data.time = (1:length(app.getTrimBeads)) / app.data.fs;
                app.applyHeader;
                app.plotGlobal;
                app.addToNew(true, [app.controls.flipDataBtn,...
                    app.controls.deselect.deselectBtn]);
                app.addToNew(true, app.controls.generateCovBtn);
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function applyHeader(app)
            % Saves important data from header and updates texts in load panel of tab 1 to
            % display this data.
            
            % Update texts.
            app.load.KATxt.String = ['K_A: ' num2str(app.load.KA) ' [pN/nm]'];
            app.load.KBTxt.String = ['K_B: ' num2str(app.load.KB) ' [pN/nm]'];
            app.load.CALATxt.String = ['CAL_A: ' num2str(app.load.CALA) ' [nm/V]'];
            app.load.CALBTxt.String = ['CAL_B: ' num2str(app.load.CALB) ' [nm/V]'];
            app.load.fsTxt.String = ['F_s: ' num2str(app.data.fs) ' [Hz]'];
            app.load.filenameTxt.String = ['File: ' app.load.name app.load.ext];
            
            % Update positioning of texts.
            drawnow;
            height = 0.04;
            app.load.CALATxt.Position = [app.load.inputBtn.Position(1) + app.load.inputBtn.Position(3) + 0.005,...
                app.load.inputBtn.Position(2) + height, app.load.CALATxt.Extent(3), height];
            app.load.CALBTxt.Position = [app.load.inputBtn.Position(1) + app.load.inputBtn.Position(3) + 0.005,...
                app.load.inputBtn.Position(2), app.load.CALBTxt.Extent(3), height];
            app.load.KATxt.Position = [app.load.CALATxt.Position(1) + app.load.CALATxt.Position(3) + 0.005,...
                app.load.inputBtn.Position(2) + height, app.load.KATxt.Extent(3), height];
            app.load.KBTxt.Position = [app.load.CALBTxt.Position(1) + app.load.CALBTxt.Position(3) + 0.005,...
                app.load.inputBtn.Position(2), app.load.KBTxt.Extent(3), height];
            app.load.fsTxt.Position = [app.load.KATxt.Position(1) + app.load.KATxt.Position(3) + 0.005,...
                app.load.inputBtn.Position(2) + height, app.load.fsTxt.Extent(3), height];
            app.load.filenameTxt.Position = [app.load.KBTxt.Position(1) + app.load.KBTxt.Position(3) + 0.005,...
                app.load.inputBtn.Position(2), app.load.filenameTxt.Extent(3), height];
            app.load.skipBtn.Position = [max([app.load.filenameTxt.Position(1) + app.load.filenameTxt.Position(3),...
                app.load.fsTxt.Position(1) + app.load.fsTxt.Position(3)]) + 0.005,...
                0.902, 0.085, 0.083];
            
            % Tell currentTask text to never overlap the skipBtn.
            app.load.taskMaxWidth = 0.995 - (app.load.skipBtn.Position(1) + app.load.skipBtn.Position(3));
            
            % Make texts visible.
            app.load.KATxt.Visible = 'on';
            app.load.KBTxt.Visible = 'on';
            app.load.CALATxt.Visible = 'on';
            app.load.CALBTxt.Visible = 'on';
            app.load.fsTxt.Visible = 'on';
            app.load.filenameTxt.Visible = 'on';
        end
        
        function plotGlobal(app)
            % Updates main axes in load panel of tab 1 with plots of beads over time.
            
            app.updateCurrentTask('Plotting bead positions.');
            
            % If beads are already plotted, xlim should be preserved.
            xliminit = [];
            if strcmp(app.load.main.bead1Plot.Visible, 'on')
                xliminit = xlim(app.load.main.mainAxes);
            end
            
            % Hide all event patches.
            for npatch = 1:length(app.load.main.patches)
                app.load.main.patches(npatch).Visible = 'off';
            end
            
            % Get beads.
            [A, B] = app.getActiveBeads;
            if app.flipped
                A = -A;
                B = -B;
            end
            
            % Update plots.
            app.changeData(app.load.main.bead1Plot, app.data.time, A);
            app.changeData(app.load.main.bead2Plot, app.data.time, B);
            app.changeData(app.controls.deselect.leftPlot, [], []);
            app.changeData(app.controls.deselect.rightPlot, [], []);
            app.load.main.mainAxes.XLabel.Visible = 'on';
            axis(app.load.main.mainAxes, 'tight');
            
            % Store outermost ylim, for plotting patches.
            app.load.main.ylim = ylim(app.load.main.mainAxes);
            
            % Update xlim, if needed.
            if ~isempty(xliminit)
                xlim(app.load.main.mainAxes, xliminit);
            end
            
            % Make visible.
            app.load.main.mainAxes.Visible = 'on';
            if strcmp(app.load.cov.covAxes.Visible, 'on')
                % If covariance axes are visible, XTickLabel and XLabel of main axes
                % should be hidden.
                app.load.main.mainAxes.XTickLabel = [];
                app.load.main.mainAxes.XTick = app.load.cov.covAxes.XTick; % Update XTicks to match.
                app.load.main.mainAxes.XLabel.Visible = 'off';
                
                % Resize axes based on yticklabels and ylabel.
                tiC = app.load.cov.covAxes.TightInset;
                posC = app.load.cov.covAxes.Position;
                tiM = app.load.main.mainAxes.TightInset;
                posM = app.load.main.mainAxes.Position;
                newWidth = min([1 - 2*posC(1) - tiC(3), 1 - 2*posM(1) - tiM(3)]);
                posC(3) = newWidth;
                posM(3) = newWidth;
                app.load.cov.covAxes.Position = posC;
                app.load.main.mainAxes.Position = posM;
            else
                % Resize axes based on yticklabels and ylabel.
                ti = app.load.main.mainAxes.TightInset;
                pos = app.load.main.mainAxes.Position;
                pos(3) = 1 - 2*pos(1) - ti(3);
                app.load.main.mainAxes.Position = pos;
            end
            
            % Update event patches, if needed.
            if ~isempty(app.load.main.patches)
                app.plotEvents;
                app.updateEventVisNColor;
            end
            
            drawnow;
        end
        
        function flipDataCallback(app, src, ~)
            % Called when user clicks flipDataBtn. Updates flipped accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.flipped = logical(src.Value);
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                if app.flipped
                    src.Value = src.Max;
                else
                    src.Value = src.Min;
                end
            end
        end
        
        function deselectCallback(app, src, ~, store, backup)
            % Called when user clicks deselectBtn. Opens the deselect panel and allows
            % user to click pairs of points on the main axes to define regions of data
            % targeted for removal.
            
            if nargin < 5
                backup = app.storeBackupTab1;
                % backup sets deselectBtn's value to what it was at the beginning of this
                % callback, which is opposite what it should be if we want to return the
                % app to its state before the button was clicked.
                backup.controls.deselect.deselectBtn.Value = ~backup.controls.deselect.deselectBtn.Value;
            end
            
            try
                switch src.Value
                    case 1
                        % If turning deselect on...
                        store = app.disable; % ...disable the app...
                        app.controls.deselect.deselectBtn.Callback = []; % ...temporarily disable deselectBtn...
                        yliminit = ylim(app.load.main.mainAxes);
                        for i = 1:length(app.controls.deselect.patches)
                            app.controls.deselect.patches(i).Visible = 'on'; % ...show all deselect patches...
                        end
                        if ~app.controls.deselect.first
                            app.controls.deselect.linePlot.Visible = 'on'; % ...as well as the red vertical line...
                        end
                        ylim(app.load.main.mainAxes, yliminit);
                        app.controls.deselect.selPatch = false(length(app.controls.deselect.patches), 1); % ...set selPatch to indicate that no patches are selected...
                        app.controls.deselect.panel.Visible = 'on';
                        app.load.main.mainAxes.XColor = app.misc.colors.activeAxes;
                        app.load.main.mainAxes.YColor = app.misc.colors.activeAxes;
                        app.load.main.mainAxes.ZColor = app.misc.colors.activeAxes;
                        app.misc.storeClickCallback = app.misc.fig.WindowButtonDownFcn;
                        app.misc.fig.WindowButtonDownFcn = {@app.deselectClick backup}; % ...set the figure button down callback to deselectClick()...
                        app.misc.storeMotionCallback = app.misc.fig.WindowButtonMotionFcn;
                        app.misc.fig.WindowButtonMotionFcn = {@app.deselectMotion backup}; % ...set the figure motion callback to deselectMotion()...
                        app.controls.deselect.deselectBtn.Callback = {@app.deselectCallback store backup}; % ...pass store to this callback...
                        app.controls.deselect.removeBtn.Callback = {@app.removeCallback store}; % ...and to removeCallback (the other way to exit this mode)...
                        app.controls.deselect.deselectBtn.Enable = 'on'; % ...enable deselectBtn...
                        app.updateCurrentTask('User is selecting sections of data to remove.'); % ...and update the current task message.
                        drawnow;
                    case 0
                        % If turning deselect off...
                        app.misc.fig.WindowButtonDownFcn = app.misc.storeClickCallback; % ...reset the figure button down callback...
                        app.misc.fig.WindowButtonMotionFcn = app.misc.storeMotionCallback; % ...as well as the figure motion callback...
                        for i = 1:length(app.controls.deselect.patches)
                            app.controls.deselect.patches(i).FaceColor = app.misc.colors.deselect; % ...reset all patches to red (i.e. not selected)...
                            app.controls.deselect.patches(i).Visible = 'off'; % ...and make them invisible...
                        end
                        app.controls.deselect.linePlot.Visible = 'off'; % ...also make red vertical line invisible...
                        app.controls.deselect.undoBtn.Enable = 'off';
                        app.controls.deselect.deselectBtn.Callback = @app.deselectCallback; % ...and reset this callback.
                        app.controls.deselect.panel.Visible = 'off';
                        app.load.main.mainAxes.XColor = 'k';
                        app.load.main.mainAxes.YColor = 'k';
                        app.load.main.mainAxes.ZColor = 'k';
                        app.restore(store);
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function undoCallback(app, ~, ~)
            % Called when user clicks undoBtn. Clears selected patches.
            
            backup = app.storeBackupTab1;
            try
                delete(app.controls.deselect.patches(app.controls.deselect.selPatch));
                app.controls.deselect.patches = app.controls.deselect.patches(~app.controls.deselect.selPatch);
                app.controls.deselect.xpos = app.controls.deselect.xpos(~app.controls.deselect.selPatch, :);
                app.controls.deselect.undoBtn.Enable = 'off';
                if size(app.controls.deselect.xpos, 1) > 0
                    app.controls.deselect.removeBtn.Enable = 'on';
                end
                app.controls.deselect.selPatch = false(length(app.controls.deselect.patches), 1);
                drawnow;
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function removeCallback(app, ~, ~, store)
            % Called when user clicks deselect.removeBtn.
            
            backup = app.storeBackupTab1;
            try
                % Display a question dialog box to the user to let them know this cannot be
                % undone.
                msg = {'The only way to restore data is to reload the original file.';'Remove anyways?'};
                promptMessage = strip(sprintf('%s\n\n%s',msg{:}));
                button = questdlg(promptMessage, 'Continue', 'Remove', 'Cancel', 'Remove');
                if strcmp(button, 'Cancel')
                    return;
                end

                app.updateCurrentTask('Removing data.');

                for i = 1:size(app.controls.deselect.xpos, 1)
                    % For each deselect patch's end coordinates...
                    seg = app.controls.deselect.xpos(i, :) * app.data.fs;
                    l = floor(seg(1)); % ...convert to left index...
                    r = ceil(seg(2)); % ...and right index.
                    if l < 1
                        if r >= 1
                            % If left index is less than 1, set it to 1 only if right index is
                            % not also less than 1.
                            l = 1;
                        else
                           continue 
                        end
                    end
                    if r > length(app.getTrimBeads)
                        if l <= length(app.getTrimBeads)
                            % If right index is greater than end, set it to end only if right
                            % index is not also greater than end.
                            r = length(app.getTrimBeads);
                        else
                            continue
                        end
                    end
                    j = find(app.controls.deselect.indices, r); % Convert from indices in trimBeads to indices in origBeads...
                    app.controls.deselect.indices(j(l):j(r)) = false; % ...and set all data points in between to be deleted.
                end

                % Save some important properties which will be deleted on reset.
                path = app.load.path;
                name = app.load.name;
                ext = app.load.ext;
                header = app.load.header;
                fs = app.data.fs;
                KA = app.load.KA;
                KB = app.load.KB;
                CALA = app.load.CALA;
                CALB = app.load.CALB;
                origData = app.data.origData;
                origBeads = app.data.origBeads;
                indices = app.controls.deselect.indices;
                simInput = app.misc.simInput;

                % Restore and reset.
                app.restore(store); % So that tab 2 elements are restored.
                app.fullReset; % Just resets tab 1.

                % Restore the above properties and run simplified version of
                % chooseTxtCallback().
                store = app.disable;
                app.load.path = path;
                app.load.name = name;
                app.load.ext = ext;
                app.load.header = header;
                app.data.fs = fs;
                app.load.KA = KA;
                app.load.KB = KB;
                app.load.CALA = CALA;
                app.load.CALB = CALB;
                app.data.origData = origData;
                app.data.origBeads = origBeads;
                app.controls.deselect.indices = indices;
                app.misc.simInput = simInput;
                app.data.time = (1:length(app.getTrimBeads)) / app.data.fs;
                
                app.updateCurrentTask('Displaying pertinent values.');
                app.applyHeader;
                
                app.plotGlobal;
                app.controls.deselect.saveBtn.Enable = 'on';
                if any(indices)
                    % Don't let user continue if they removed all data.
                    app.addToNew(true, [app.controls.flipDataBtn,...
                        app.controls.deselect.deselectBtn]);
                    app.addToNew(true, app.controls.generateCovBtn);
                end
                
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function saveCallback(app, ~, ~)
            % Called when user clicks saveBtn. Creates a new .txt file with header
            % matching the current header and data replaced by trimData.
            
            backup = app.storeBackupTab1;
            try
                % Let user specify name for new file.
                tempName = [app.load.name ' CUT' app.load.ext];
                [name, path] = uiputfile([app.load.path '/' tempName]);
                if name == 0
                    return
                end
                [path, name] = fileparts([path name]);
                name = [name '.txt'];

                app.updateCurrentTask(sprintf('Saving data to %s.', fullfile(path, name)));

                app.controls.deselect.saveBtn.Enable = 'off';
                drawnow;

                % Write header and data to file.
                try
                    sHeader = app.load.txtOutputFcn(...
                        app.data.fs, app.load.KA,...
                        app.load.KB, app.load.CALA,...
                        app.load.CALB, app.data.origData,...
                        app.load.header);
                catch ME
                    % Display any error to user and exit.
                    app.unexpectedError(ME, sprintf('There was an error while writing your data:\n\n%s', ME.message));
                    app.controls.deselect.saveBtn.Enable = 'on';
                    app.updateCurrentTask('User is selecting sections of data to remove.');
                    return
                end
                fileID = fopen(fullfile(path, name), 'w');
                if fileID ~= -1
                    fprintf(fileID, '%s', sHeader);
                    aData = table2array(app.getTrimData);
                    numCols = size(aData, 2);
                    format = [repmat('%f\t', 1, numCols-1) '%f\n'];
                    fprintf(fileID, format, aData');
                    fclose(fileID);
                end

                % Update program so that the newly made file replaces the current file.
                [app.load.path, app.load.name, app.load.ext] = fileparts(fullfile(path, name));
                app.load.filenameTxt.String = ['File: ' app.load.name app.load.ext];
                app.load.filenameTxt.Position(3) = app.load.filenameTxt.Extent(3);
                app.data.origData = app.getTrimData;
                app.data.origBeads = app.getTrimBeads;
                app.controls.deselect.indices = true(size(app.data.origBeads,1), 1);

                app.updateCurrentTask('User is selecting sections of data to remove.');
            catch ME
                try
                    fclose(fileID);
                catch
                end
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME, sprintf('An unexpected error occurred while saving:\n\n%s', ME.message));
                app.restoreBackupTab1(backup);
            end
        end
        
        function deselectClick(app, ~, ~, backup)
            % Button down callback for the figure when user is in deselect mode. Allows
            % user to click on the main axes to place deselect patches which target data
            % to be removed.
            
            try
                if ~isequal(app.misc.tab1, app.misc.tabgroup.SelectedTab)
                    return
                end
                limx = app.load.main.mainAxes.XLim;
                limy = app.load.main.mainAxes.YLim;
                ex = 0.0025*(limx(2)-limx(1));
                [in, p1] = app.inAxes(app.load.main.mainAxes); % [whether in axes, mouse position in axes]
                if abs(p1(1) - limx(1)) < ex && p1(2) > limy(1) && p1(2) < limy(2)
                    % If mouse X coordinate is close to lower xlim and mouse Y coordinate is
                    % within ylim...
                    x = app.controls.deselect.leftPlot.XData(1); % ...set x to leftPlot XData.
                elseif abs(p1(1) - limx(2)) < ex && p1(2) > limy(1) && p1(2) < limy(2)
                    % If mouse X coordinate is close to upper xlim and mouse Y coordinate is
                    % within ylim...
                    x = app.controls.deselect.rightPlot.XData(1); % ...set x to rightPlot XData.
                elseif in
                    % Otherwise, if mouse is in the main axes...
                    x = p1(1); % ...set x to mouse X coordinate.
                else
                    return
                end
                limy = app.load.main.ylim;
                % Clicks come in pairs, with each click defining one edge of the patch.
                if app.controls.deselect.first
                    % If this click marks the first click in a pair...
                    if ~isempty(app.controls.deselect.xpos)
                        for i = 1:length(app.controls.deselect.patches)
                            app.controls.deselect.patches(i).FaceColor = app.misc.colors.deselect;
                        end
                        app.controls.deselect.undoBtn.Enable = 'off';
                        inPrev = app.controls.deselect.xpos(:,1) <= x & x <= app.controls.deselect.xpos(:,2);
                        if any(inPrev)
                            % ...it can fall within a previous patch...
                            if app.controls.deselect.selPatch == inPrev
                                % If this previous patch was already selected, unselect it.
                                app.controls.deselect.patches(inPrev).FaceColor = app.misc.colors.deselect; % (make it red)
                                app.controls.deselect.selPatch = false(length(app.controls.deselect.patches), 1); % (update selPatch)
                                app.controls.deselect.removeBtn.Enable = 'on'; % (let user remove data)
                            else
                                % If this previous patch was not already selected, select it.
                                app.controls.deselect.patches(inPrev).FaceColor = app.misc.colors.deselectSelected; % (make it blue)
                                app.controls.deselect.selPatch = inPrev; % (update selPatch)
                                app.controls.deselect.undoBtn.Enable = 'on'; % (let user undo the patch)
                                app.controls.deselect.removeBtn.Enable = 'off';
                            end
                            return
                        end
                    end
                    % ...or it can mark the beginning of a new patch.
                    app.controls.deselect.removeBtn.Enable = 'off';
                    app.changeData(app.controls.deselect.linePlot, [x x], limy); % (indicate this with the red vertical line)
                    app.controls.deselect.xpos(end+1,1) = x; % (save the location to xpos)
                    app.controls.deselect.selPatch = false(length(app.controls.deselect.patches), 1); % (update selPatch)
                else
                    % If this click marks the second click in a pair...
                    app.controls.deselect.xpos(end,2) = x; % ...save the location to xpos.

                    % Determine if the new patch overlaps with any old patches.
                    all = app.controls.deselect.xpos;
                    allbut = all(1:end-1,:); % The old patches.
                    sel = sort(all(end,:)); % The new patch.
                    dellist = false(size(allbut, 1), 1);
                    for i = 1:size(allbut, 1)
                        % For each old patch...
                        prev = allbut(i, :);
                        if sel(1) >= prev(1) && sel(1) <= prev(2)
                            % ...if the left boundary of the new patch falls within the old
                            % patch...
                            sel(1) = prev(1); % ...set the left boundary of the new patch to the left boundary of the old patch...
                            dellist(i) = true; % ...and delete the old patch.
                        end
                        if sel(2) >= prev(1) && sel(2) <= prev(2)
                            % ...if the right boundary of the new patch falls within the old
                            % patch...
                            sel(2) = prev(2); % ...set the right boundary of the new patch to the right boundary of the old patch...
                            dellist(i) = true; % ...and delete the old patch.
                        end
                        if sel(1) <= prev(1) && sel(2) >= prev(2)
                            % ...if the new patch encompasses the old patch...
                            dellist(i) = true; % ...delete the old patch.
                        end
                    end
                    all = [allbut(~dellist,:); sel];
                    app.controls.deselect.xpos = all;
                    oldPatches = app.controls.deselect.patches;

                    % Plot the new patch.
                    yliminit = ylim(app.load.main.mainAxes);
                    limy(1) = limy(1) - 0.1*(limy(2) - limy(1));
                    limy(2) = limy(2) + 0.1*(limy(2) - limy(1));
                    app.controls.deselect.patches(end+1) = patch(app.load.main.mainAxes,...
                        'XData', [sel(1) sel(2) sel(2) sel(1)],...
                        'YData', [limy(1) limy(1) limy(2) limy(2)],...
                        'FaceColor', app.misc.colors.deselect,...
                        'FaceAlpha', 0.3,...
                        'HitTest', 'off');
                    ylim(app.load.main.mainAxes, yliminit);
                    
                    % Remove the old patch(es) and add the new patch to list of patches.
                    delete(oldPatches(dellist));
                    app.controls.deselect.patches = [oldPatches(~dellist) app.controls.deselect.patches(end)];
                    app.controls.deselect.selPatch = false(length(app.controls.deselect.patches), 1);
                    app.controls.deselect.linePlot.Visible = 'off';
                    app.controls.deselect.removeBtn.Enable = 'on';
                end
                app.controls.deselect.first = ~app.controls.deselect.first;
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowButtonDownFcn = [];
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function deselectMotion(app, ~, ~, backup)
            % Motion callback for the figure when user is in deselect mode. Highlights
            % lower and upper xlim of main axes with black lines when user hovers over
            % them.
            
            thisMotionCallback = app.misc.fig.WindowButtonMotionFcn;
            try
                if ~isequal(app.misc.tab1, app.misc.tabgroup.SelectedTab)
                    return
                end
                limx = app.load.main.mainAxes.XLim;
                limy = app.load.main.mainAxes.YLim;
                ex = 0.0025*(limx(2)-limx(1));
                [in, p1] = app.inAxes(app.load.main.mainAxes, ex); % [whether in axes, mouse position in axes]
                if in
                    % If mouse is in the axes, set pointer to crosshair.
                    app.misc.fig.Pointer = 'crosshair';
                else
                    % Otherwise, reset pointer to arrow.
                    app.misc.fig.Pointer = 'arrow';
                end
                if abs(p1(1) - limx(1)) < ex && p1(2) > limy(1) && p1(2) < limy(2)
                    % If mouse X coordinate is close to lower xlim and mouse Y coordinate is
                    % within ylim...
                    app.changeData(app.controls.deselect.leftPlot, [limx(1) limx(1)], limy);
                    app.controls.deselect.leftPlot.Color = 'k'; % ...show black line at lower xlim.
                elseif abs(p1(1) - limx(2)) < ex && p1(2) > limy(1) && p1(2) < limy(2)
                    % If mouse X coordinate is close to upper xlim and mouse Y coordinate is
                    % within ylim...
                    app.changeData(app.controls.deselect.rightPlot, [limx(2) limx(2)], limy);
                    app.controls.deselect.rightPlot.Color = 'k'; % ...show black line at upper xlim.
                else
                    app.controls.deselect.rightPlot.Color = 'none';
                    app.controls.deselect.leftPlot.Color = 'none';
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                if isequal(app.misc.fig.WindowButtonMotionFcn, thisMotionCallback)
                    app.misc.fig.WindowButtonMotionFcn = [];
                    app.unexpectedError(ME);
                    app.restoreBackupTab1(backup);
                end
            end
        end
        
        function generateCovCallback(app, ~, ~)
            % Called when user clicks generateCovBtn. Calls autoGenerateCov() and enables
            % controls panel of tab 1.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.autoGenerateCov;
                app.addToNew(false, app.controls.generateCovBtn);
                app.addToNew(true, [app.controls.windows.covwindowTxt,...
                    app.controls.windows.covwindowEdit,...
                    app.controls.windows.covsmoothTxt,...
                    app.controls.windows.covsmoothEdit,...
                    app.controls.windows.interactBtn,...
                    app.controls.windows.resetBtn,...
                    app.controls.cutoffs.minDurTxt,...
                    app.controls.cutoffs.minDurEdit,...
                    app.controls.cutoffs.minSepTxt,...
                    app.controls.cutoffs.minSepEdit,...
                    app.controls.cutoffs.resetBtn,...
                    app.controls.useChangepointBtn,...
                    app.controls.peakToPeakBtn,...
                    app.controls.thresholdBtn,...
                    app.controls.automaticBtn,...
                    app.controls.manualBtn,...
                    app.controls.setManualBtn,...
                    app.controls.peak1Txt,...
                    app.controls.peak1Edit,...
                    app.controls.peak2Txt,...
                    app.controls.peak2Edit,...
                    app.controls.minTxt,...
                    app.controls.minEdit,...
                    app.controls.beadABtn,...
                    app.controls.beadBBtn,...
                    app.controls.bothBeadsBtn]);
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function autoGenerateCov(app)
            % Calculates the covariance, updates the covariance and covariance histogram
            % plots, updates the automatic peaks/min based on the histogram plot, and
            % updates the preliminary event patches, if needed.
            
            app.updateCurrentTask('Calculating the covariance.');
            covwin = app.controls.windows.covwindow;
            covsm = app.controls.windows.covsmooth;
            [A, B] = app.getActiveBeads;
            app.load.cov.cov = app.calculateCov(A, B, covwin, covsm, app.misc.signal);
            
            app.updateCurrentTask('Plotting the covariance.');
            app.plotCov(app.load.cov.cov);
            
            app.updateCurrentTask('Plotting the covariance histogram.');
            app.plotCovHist(app.load.cov.cov);
            
            app.updateCurrentTask('Determining the peaks and minimum.');
            edges = app.analyze.covHist.covHistPlot.BinEdges;
            values = app.analyze.covHist.covHistPlot.Values;
            app.autop1 = []; app.autop2 = []; app.automin = [];
            [app.autop1, app.autop2, app.automin] = app.calculatePeaksMin(edges, values);
            if ~isempty(app.misc.listenerError)
                rethrow(app.misc.listenerError)
            end
            
            if ~isempty(app.load.cov.patches)
                % Only update if already showing.
                cov = app.load.cov.cov;
                switch app.useman
                    case false
                        p1 = app.autop1;
                        p2 = app.autop2;
                        m = app.automin;
                    case true
                        p1 = app.manp1;
                        p2 = app.manp2;
                        m = app.manmin;
                end
                switch app.usemin
                    case false
                        prelimEvents = app.findPrelimEvents(cov, p1, p2);
                    case true
                        prelimEvents = app.findPrelimEventsByMin(cov, m);
                end
                app.plotPrelimEvents(prelimEvents);
            end
            
            % If manual peaks/min have not been set by the user, set them equal to their
            % automatic counterparts.
            if isempty(app.manmin)
                app.manmin = app.automin;
            end
            if isempty(app.manp1)
                app.manp1 = app.autop1;
            end
            if isempty(app.manp2)
                app.manp2 = app.autop2;
            end
            if ~isempty(app.misc.listenerError)
                rethrow(app.misc.listenerError)
            end
            
            drawnow;
        end
        
        function covwindowCallback(app, src, ~)
            % Called when user changes the value of covwindowEdit. Makes sure input is
            % valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                if ~isfinite(d) || ~isreal(d) || d <= 0
                    % If input is infinite, NaN, or not positive, reset edit text string to
                    % current value and return.
                    uiwait(errordlg('The averaging window must be a positive integer.', 'Error Dialog', 'modal'))
                    src.String = int2str(app.controls.windows.covwindow);
                    app.restore(store);
                    return
                end

                if mod(d,1) ~= 0
                    % If input is a decimal, round down to nearest integer.
                    warndlg('The averaging window must be an integer. Rounding input down.')
                    d = floor(d);
                end
                % Input is a positive integer here.

                % If new value does not match current value, covariance will be different.
                change = app.controls.windows.covwindow ~= d;

                % Update edit text and current value.
                src.String = int2str(d);
                app.controls.windows.covwindow = d;

                if change
                    % If covariance will be different, recalculate it and update plots.
                    app.autoGenerateCov;
                end
                
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.String = int2str(app.controls.windows.covwindow);
            end
        end
        
        function covsmoothCallback(app, src, ~)
            % Called when user changes the value of covsmoothEdit. Makes sure input is
            % valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                if ~isfinite(d) || ~isreal(d) || d <= 0
                    % If input is infinite, NaN, or not positive, reset edit text string to
                    % current value and return.
                    uiwait(errordlg('The smoothing window must be a positive integer.', 'Error Dialog', 'modal'))
                    src.String = int2str(app.controls.windows.covsmooth);
                    app.restore(store);
                    return
                end

                if mod(d,1) ~= 0
                    % If input is a decimal, round down to nearest integer.
                    warndlg('The smoothing window must be an integer. Rounding input down.')
                    d = floor(d);
                end

                % Input is a positive integer here.
                d = 2*round((d+1)/2)-1;
                % Input is a positive odd integer here.

                if d < 3
                    warndlg('The smoothing window must be greater than the order of the filter, which is 2. Setting it to 3.')
                    d = 3;
                end

                % If new value does not match current value, covariance will be different.
                change = app.controls.windows.covsmooth ~= d;

                % Update edit text and current value.
                src.String = int2str(d);
                app.controls.windows.covsmooth = d;

                if change
                    % If covariance will be different, recalculate it and update plots.
                    app.autoGenerateCov;
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.String = int2str(app.controls.windows.covsmooth);
            end
        end
        
        function windowsCallback(app, ~, ~)
            % Called when user clicks interactBtn. Opens the interactive filter window
            % chooser.
            
            backup = app.storeBackupTab1;
            try
                % Temporarily disable click callback so user doesn't click before store is
                % passed to the callback.
                app.misc.storeClickCallback = app.misc.fig.WindowButtonDownFcn;
                app.misc.fig.WindowButtonDownFcn = [];

                store = app.disable;
                drawnow;
                app.controls.windows.panel.Visible = 'on'; % Show interactive axes.
                app.changeData(app.controls.windows.pt, app.controls.windows.covwindow, app.controls.windows.covsmooth); % Place a point at current covwindow and covsmooth.
                app.useman = false; % Set program to automatically determine peaks/min.
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                app.misc.storeMotionCallback = app.misc.fig.WindowButtonMotionFcn;
                app.misc.fig.WindowButtonMotionFcn = {@app.windowsMotion backup}; % Set the figure motion callback to windowsMotion().
                app.misc.fig.WindowButtonDownFcn = {@app.windowsClick store backup}; % Set the figure button down callback to windowsClick(), and pass store.
                app.updateCurrentTask('User is choosing windows interactively.'); % Update the current task message.
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function windowsMotion(app, ~, ~, backup)
            % Motion callback for the figure when interactive filter window chooser is
            % open.
            
            thisMotionCallback = app.misc.fig.WindowButtonMotionFcn;
            try
                [in, p1] = app.inAxes(app.controls.windows.windowsAxes); % [whether in axes, mouse position in axes]
                if in
                    % If in the axes, set point's position to mouse X and mouse Y...
                    app.changeData(app.controls.windows.pt, p1(1), p1(2));
                else
                    % Otherwise...
                    if app.controls.windows.pt.XData == app.controls.windows.covwindow && app.controls.windows.pt.YData == app.controls.windows.covsmooth
                        % ...if point's position is up to date, reset the click callback and return.
                        return
                    end
                    % ...if point's position is outdated, set it to current covwindow and
                    % covsmooth...
                    app.changeData(app.controls.windows.pt, app.controls.windows.covwindow, app.controls.windows.covsmooth);
                end
                % ...and update covariance and covariance histogram plots.
                covwindow = round(app.controls.windows.pt.XData);
                covsmooth = 2*round((app.controls.windows.pt.YData+1)/2)-1;
                [A, B] = app.getActiveBeads;
                cov = app.calculateCov(A, B, covwindow, covsmooth, app.misc.signal);
                app.plotCov(cov);
                app.plotCovHist(cov);
                edges = app.analyze.covHist.covHistPlot.BinEdges;
                values = app.analyze.covHist.covHistPlot.Values;
                [p1, p2, m] = app.calculatePeaksMin(edges, values);
                if app.usemin
                    p1 = [];
                    p2 = [];
                else
                    m = [];
                end
                app.updatePeaksMinPlots(p1, p2, m);
                
                % Update edit texts.
                app.controls.windows.covwindowEdit.String = int2str(covwindow);
                app.controls.windows.covsmoothEdit.String = int2str(covsmooth);
                app.controls.peak1Edit.String = num2str(p1);
                app.controls.peak2Edit.String = num2str(p2);
                app.controls.minEdit.String = num2str(m);
            catch ME
                if ~isvalid(app)
                    return
                end
                if isequal(app.misc.fig.WindowButtonMotionFcn, thisMotionCallback)
                    app.misc.fig.WindowButtonMotionFcn = [];
                    app.unexpectedError(ME);
                    app.restoreBackupTab1(backup);
                end
            end
        end
        
        function windowsClick(app, ~, ~, store, backup)
            % Button down callback for the figure when interactive filter window chooser
            % is open. When user clicks, regardless of mouse position, covwindow and
            % covsmooth are updated based on point's position and plots are updated
            % accordingly.
            
            try
                app.controls.windows.panel.Visible = 'off'; % Hide interactive axes.
                app.misc.fig.WindowButtonMotionFcn = app.misc.storeMotionCallback; % Reset the figure motion callback.
                app.misc.fig.WindowButtonDownFcn = app.misc.storeClickCallback; % Reset the figure button down callback.
                drawnow;
                % Update values.
                app.controls.windows.covwindow = round(app.controls.windows.pt.XData);
                app.controls.windows.covsmooth = 2*round((app.controls.windows.pt.YData+1)/2)-1;

                % Update edit texts.
                app.controls.windows.covwindowEdit.String = int2str(app.controls.windows.covwindow);
                app.controls.windows.covsmoothEdit.String = int2str(app.controls.windows.covsmooth);

                % Update plots.
                app.autoGenerateCov;
                
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowButtonDownFcn = [];
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function defwindowsCallback(app, ~, ~)
            % Called when user clicks windows.resetBtn. Restores covwindow and covsmooth
            % to default values and updates plots, if needed.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % If either covwindow or covsmooth is different from its old value, covariance
                % will be different.
                change = app.controls.windows.covwindow ~= app.controls.windows.defcovwindow ||...
                    app.controls.windows.covsmooth ~= app.controls.windows.defcovsmooth;

                % Update values.
                app.controls.windows.covwindow = app.controls.windows.defcovwindow;
                app.controls.windows.covsmooth = app.controls.windows.defcovsmooth;

                % Update edit texts.
                app.controls.windows.covwindowEdit.String = int2str(app.controls.windows.covwindow);
                app.controls.windows.covsmoothEdit.String = int2str(app.controls.windows.covsmooth);

                if change
                    app.autoGenerateCov;
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function plotCov(app, cov)
            % Updates covariance axes in load panel of tab 1 with plot of covariance over
            % time.
            
            % If this is the first time covariance is being plotted...
            init = strcmp(app.load.cov.covPlot.Visible, 'off');
            
            % Clear plots of peaks/min and preliminary event patches.
            app.changeData(app.load.cov.peak1Plot, [], []);
            app.changeData(app.load.cov.peak2Plot, [], []);
            app.changeData(app.load.cov.minPlot, [], []);
            for npatch = 1:length(app.load.cov.patches)
                app.load.cov.patches(npatch).Visible = 'off';
            end
            
            % Plot covariance.
            app.changeData(app.load.cov.covPlot, app.data.time, cov);
            axis(app.load.cov.covAxes, 'tight');
            
            if init
                % ...then plot should be made visible and XTickLabel and XLabel of main
                % axes should be hidden.
                app.load.cov.covAxes.Visible = 'on';
                app.load.main.mainAxes.XTickLabel = [];
                app.load.main.mainAxes.XTick = app.load.cov.covAxes.XTick;
                app.load.main.mainAxes.XLabel.Visible = 'off';
                linkaxes([app.load.main.mainAxes, app.load.cov.covAxes], 'x');
            end
            
            % Resize axes based on yticklabels and ylabel.
            tiC = app.load.cov.covAxes.TightInset;
            posC = app.load.cov.covAxes.Position;
            tiM = app.load.main.mainAxes.TightInset;
            posM = app.load.main.mainAxes.Position;
            newWidth = min([1 - 2*posC(1) - tiC(3), 1 - 2*posM(1) - tiM(3)]);
            posC(3) = newWidth;
            posM(3) = newWidth;
            app.load.cov.covAxes.Position = posC;
            app.load.main.mainAxes.Position = posM;
        end
        
        function plotCovHist(app, cov)
            % Updates covariance histogram axes in analyze panel of tab 1 with plot of
            % covariance histogram.
            
            % If this is the first time covariance histogram is being plotted...
            init = strcmp(app.analyze.covHist.covHistPlot.Visible, 'off');
            
            % Clear plots of peaks/min.
            app.changeData(app.analyze.covHist.peak1Plot, [], []);
            app.changeData(app.analyze.covHist.peak2Plot, [], []);
            app.changeData(app.analyze.covHist.minPlot, [], []);
            
            % Plot covariance histogram.
            app.analyze.covHist.covHistPlot.Data = cov;
            app.analyze.covHist.covHistPlot.NumBins = 100;
            axis(app.analyze.covHist.covHistAxes, 'auto xy');
            
            if init
                % ...then plot should be made visible.
                app.analyze.covHist.covHistPlot.Visible = 'on';
                app.analyze.covHist.covHistAxes.Visible = 'on';
            end
            
            app.analyze.covHist.covHistAxes.YAxis.Exponent = 0;
            
            % Resize axes based on yticklabels and ylabel.
            ti = app.analyze.covHist.covHistAxes.TightInset;
            pos = app.analyze.covHist.covHistAxes.Position;
            pos(3) = 0.33 - pos(1) - ti(3);
            app.analyze.covHist.covHistAxes.Position = pos;
        end
        
        function updatePeaksMinPlots(app, p1, p2, m)
            % Updates covariance axes in load panel of tab 1 and covariance histogram axes
            % in analyze panel of tab 1 with plots of peak 1 and peak 2 or minimum,
            % depending on usemin.
            
            if isempty(p1)
                % A p1 value of [] indicates that the peak 1 plots should be hidden.
                app.changeData(app.analyze.covHist.peak2Plot, [], []);
                app.changeData(app.load.cov.peak2Plot, [], []);
            else
                if strcmp(app.analyze.covHist.covHistAxes.Visible, 'on')
                    % As long as the covariance histogram plot is showing, update the peak
                    % 1 line.
                    limy = ylim(app.analyze.covHist.covHistAxes);
                    app.changeData(app.analyze.covHist.peak1Plot, [p1 p1], limy);
                end
                if strcmp(app.load.cov.covAxes.Visible, 'on')
                    % As long as the covariance plot is showing, update the peak 1 line.
                    limx = xlim(app.load.cov.covAxes);
                    app.changeData(app.load.cov.peak1Plot, limx, [p1 p1]);
                end
            end
            
            if isempty(p2)
                % A p2 value of [] indicates that the peak 2 plots should be hidden.
                app.changeData(app.analyze.covHist.peak1Plot, [], []);
                app.changeData(app.load.cov.peak1Plot, [], []);
            else
                if strcmp(app.analyze.covHist.covHistAxes.Visible, 'on')
                    % As long as the covariance histogram plot is showing, update the peak
                    % 2 line.
                    limy = ylim(app.analyze.covHist.covHistAxes);
                    app.changeData(app.analyze.covHist.peak2Plot, [p2 p2], limy);
                end
                if strcmp(app.load.cov.covAxes.Visible, 'on')
                    % As long as the covariance plot is showing, update the peak 2 line.
                    limx = xlim(app.load.cov.covAxes);
                    app.changeData(app.load.cov.peak2Plot, limx, [p2 p2]);
                end
            end
            
            if isempty(m)
                % A m value of [] indicates that the minimum plots should be hidden.
                app.changeData(app.analyze.covHist.minPlot, [], []);
                app.changeData(app.load.cov.minPlot, [], []);
            else
                if strcmp(app.analyze.covHist.covHistAxes.Visible, 'on')
                    % As long as the covariance histogram plot is showing, update the
                    % minimum line.
                    limy = ylim(app.analyze.covHist.covHistAxes);
                    app.changeData(app.analyze.covHist.minPlot, [m m], limy);
                end
                if strcmp(app.load.cov.covAxes.Visible, 'on')
                    % As long as the covariance plot is showing, update the minimum line.
                    limx = xlim(app.load.cov.covAxes);
                    app.changeData(app.load.cov.minPlot, limx, [m m]);
                end
            end
        end
        
        function detectionMethodCallback(app, src, eventData)
            % Called when user clicks either radio button in eventDetectMethodBtnGrp.
            % Updates usemin accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.usemin = logical(app.controls.thresholdBtn.Value);
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.SelectedObject = eventData.OldValue;
            end
        end
        
        function manvautoCallback(app, src, eventData)
            % Called when user clicks either radio button in peakDetectMethodBtnGrp.
            % Updates useman accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.useman = logical(app.controls.manualBtn.Value);
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.SelectedObject = eventData.OldValue;
            end
        end
        
        function setManualCallback(app, ~, ~)
            % Called when user clicks manualBtn. Allows user to then click three points on
            % the covariance histogram axes which will mark peak 1, peak 2, and minimum.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Clear plots of peaks/min on covariance histogram axes and covariance axes.
                app.changeData(app.analyze.covHist.peak1Plot, [], []);
                app.changeData(app.analyze.covHist.peak2Plot, [], []);
                app.changeData(app.analyze.covHist.minPlot, [], []);
                app.changeData(app.load.cov.peak1Plot, [], []);
                app.changeData(app.load.cov.peak2Plot, [], []);
                app.changeData(app.load.cov.minPlot, [], []);

                % Clear edit texts.
                app.controls.peak1Edit.String = '';
                app.controls.peak2Edit.String = '';
                app.controls.minEdit.String = '';

                app.misc.storeMotionCallback = app.misc.fig.WindowButtonMotionFcn;
                app.misc.fig.WindowButtonMotionFcn = {@app.manualMotion backup}; % Set the figure motion callback to manualMotion().
                app.misc.storeClickCallback = app.misc.fig.WindowButtonDownFcn;
                app.misc.fig.WindowButtonDownFcn = {@app.manualClick store backup []}; % Set the figure button down callback to manualClick(), and pass store.
                app.misc.fig.WindowKeyPressFcn = {@app.manualEscape store backup}; % Let user escape by pressing escape key.
                app.updateCurrentTask('User is setting the locations of the peaks/min.'); % Update the current task message.

                app.analyze.covHist.covHistAxes.XColor = app.misc.colors.activeAxes;
                app.analyze.covHist.covHistAxes.YColor = app.misc.colors.activeAxes;
                app.analyze.covHist.covHistAxes.ZColor = app.misc.colors.activeAxes;
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
                
        function manualMotion(app, ~, ~, backup)
            % Motion callback for the figure when user is setting the manual peaks/min by
            % clicking the covariance histogram plot.
            
            thisMotionCallback = app.misc.fig.WindowButtonMotionFcn;
            try
                [in, p1] = app.inAxes(app.analyze.covHist.covHistAxes); % [whether in axes, mouse position in axes]
                if isequal(app.misc.tab1, app.misc.tabgroup.SelectedTab)
                    if in
                        % Plot a red vertical line at mouse X position.
                        limy = ylim(app.analyze.covHist.covHistAxes);
                        app.changeData(app.analyze.covHist.minPlot, [p1(1) p1(1)], limy);
                        limx = xlim(app.load.cov.covAxes);
                        app.changeData(app.load.cov.minPlot, limx, [p1(1) p1(1)]);
                    else
                        app.changeData(app.analyze.covHist.minPlot, [], []);
                        app.changeData(app.load.cov.minPlot, [], []);
                    end
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                if isequal(app.misc.fig.WindowButtonMotionFcn, thisMotionCallback)
                    app.misc.fig.WindowButtonMotionFcn = [];
                    app.unexpectedError(ME);
                    app.restoreBackupTab1(backup);
                end
            end
        end
        
        function manualClick(app, ~, ~, store, backup, out)
            % Button down callback for the figure when user is setting the manual
            % peaks/min by clicking the covariance histogram plot. out contains the X
            % coordinates of the user's previous clicks.
            
            try
                if ~isequal(app.misc.tab1, app.misc.tabgroup.SelectedTab)
                    return
                end

                fcn = app.misc.fig.WindowButtonDownFcn;
                app.misc.fig.WindowButtonDownFcn = []; % Temporarily disable button down callback.

                [in, p1] = app.inAxes(app.analyze.covHist.covHistAxes); % [whether in axes, mouse position in axes]
                if in
                    out(end+1) = p1(1); % Add mouse X coordinate to out array.
                    limy = ylim(app.analyze.covHist.covHistAxes);
                    limx = xlim(app.load.cov.covAxes);
                    switch length(out)
                        case 1
                            % First click.
                            app.changeData(app.analyze.covHist.peak1Plot, [p1(1) p1(1)], limy); % Plot a red vertical line at mouse X position.
                            app.changeData(app.load.cov.peak1Plot, limx, [p1(1) p1(1)]);
                            app.misc.fig.WindowButtonDownFcn = {@app.manualClick store backup out}; % Pass out to next click's callback.
                        case 2
                            % Second click.
                            app.changeData(app.analyze.covHist.peak2Plot, [p1(1) p1(1)], limy); % Plot a red vertical line at mouse X position.
                            app.changeData(app.load.cov.peak2Plot, limx, [p1(1) p1(1)]);
                            app.misc.fig.WindowButtonDownFcn = {@app.manualClick store backup out}; % Pass out to next click's callback.
                        case 3
                            % Third click.
                            app.misc.fig.WindowButtonMotionFcn = app.misc.storeMotionCallback; % Reset the figure motion callback.
                            app.misc.fig.WindowButtonDownFcn = app.misc.storeClickCallback; % Reset the figure button down callback.
                            app.misc.fig.WindowKeyPressFcn = @app.nullKeyPressCallback; % Reset the figure key press callback.

                            % Sort three X coordinates to determine peak 1, peak 2, and
                            % minimum, and set program to use manually determined peaks/mins.
                            app.manp1 = min(out);
                            app.manp2 = max(out);
                            app.manmin = median(out);
                            app.useman = true;
                            if ~isempty(app.misc.listenerError)
                                rethrow(app.misc.listenerError)
                            end
                            
                            app.analyze.covHist.covHistAxes.XColor = 'k';
                            app.analyze.covHist.covHistAxes.YColor = 'k';
                            app.analyze.covHist.covHistAxes.ZColor = 'k';

                            app.restore(store);
                    end
                    return
                end

                % If user clicked outside of axes, reset button down callback to what it was
                % previously.
                app.misc.fig.WindowButtonDownFcn = fcn;
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowButtonDownFcn = [];
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function manualEscape(app, ~, eventData, store, backup)
            % When set manual mode is active, this is called whenever user presses any
            % key. If user presses the escape key, program will exit set manual mode and
            % return to state prior to entering the mode.
            
            try
                if strcmp(eventData.Key, 'escape')
                    app.misc.fig.WindowButtonMotionFcn = app.misc.storeMotionCallback; % Reset the figure motion callback.
                    app.misc.fig.WindowButtonDownFcn = app.misc.storeClickCallback; % Reset the figure button down callback.
                    app.misc.fig.WindowKeyPressFcn = @app.nullKeyPressCallback; % Reset the figure key press callback.

                    % Leave peaks/min as they were, but tell program to redraw.
                    app.useman = ~app.useman;
                    app.useman = ~app.useman;
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end

                    % Reset axes colors.
                    app.analyze.covHist.covHistAxes.XColor = 'k';
                    app.analyze.covHist.covHistAxes.YColor = 'k';
                    app.analyze.covHist.covHistAxes.ZColor = 'k';

                    app.restore(store);
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowKeyPressFcn = [];
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function nullKeyPressCallback(~, ~, ~)
            % Do nothing. This prevents figure from losing focus on key presses.
        end
        
        function peak1Callback(app, src, ~)
            % Called when user changes the value of peak1Edit. Makes sure input is valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                % Get current value.
                switch app.useman
                    case false
                        p1 = app.autop1;
                    case true
                        p1 = app.manp1;
                end

                if ~isfinite(d) || ~isreal(d)
                    % If input is infinite or NaN, reset edit text string to current value and
                    % return.
                    uiwait(errordlg('Value for peak 1 should be numeric and finite.', 'Error Dialog', 'modal'))
                    src.String = num2str(p1);
                    app.restore(store);
                    return
                end

                if p1 ~= d
                    % If new value does not match current value, update current value and set
                    % program to use manually determined peaks/min.
                    app.manp1 = d;
                    app.useman = true;
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                switch app.useman
                    case false
                        p1 = app.autop1;
                    case true
                        p1 = app.manp1;
                end
                src.String = num2str(p1);
            end
        end
        
        function peak2Callback(app, src, ~)
            % Called when user changes the value of peak2Edit. Makes sure input is valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                % Get current value.
                switch app.useman
                    case false
                        p2 = app.autop2;
                    case true
                        p2 = app.manp2;
                end

                if ~isfinite(d) || ~isreal(d)
                    % If input is infinite or NaN, reset edit text string to current value and
                    % return.
                    uiwait(errordlg('Value for peak 2 should be numeric and finite.', 'Error Dialog', 'modal'))
                    src.String = num2str(p2);
                    app.restore(store);
                    return
                end

                if p2 ~= d
                    % If new value does not match current value, update current value and set
                    % program to use manually determined peaks/min.
                    app.manp2 = d;
                    app.useman = true;
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                switch app.useman
                    case false
                        p2 = app.autop2;
                    case true
                        p2 = app.manp2;
                end
                src.String = num2str(p2);
            end
        end
        
        function minCallback(app, src, ~)
            % Called when user changes the value of minEdit. Makes sure input is valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                % Get current value.
                switch app.useman
                    case false
                        m = app.automin;
                    case true
                        m = app.manmin;
                end

                if ~isfinite(d) || ~isreal(d)
                    % If input is infinite or NaN, reset edit text string to current value and
                    % return.
                    uiwait(errordlg('Value for the minimum threshold should be numeric and finite.', 'Error Dialog', 'modal'))
                    src.String = num2str(m);
                    app.restore(store);
                    return
                end

                if m ~= d
                    % If new value does not match current value, update current value and set
                    % program to use manually determined peaks/min.
                    app.manmin = d;
                    app.useman = true;
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end
                end
                
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                switch app.useman
                    case false
                        m = app.automin;
                    case true
                        m = app.manmin;
                end
                src.String = num2str(m);
            end
        end
        
        function whichBeadCallback(app, src, eventData)
            % Called when user clicks any of the radio buttons in the whichBeadBtnGrp.
            % Checks if the events need to be updated.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.checkOutdatedEvents;
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.SelectedObject = eventData.OldValue;
            end
        end
        
        function useChangepointCallback(app, src, ~)
            % Called when user clicks useChangepointBtn. Updates useChangepoint
            % accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.checkOutdatedEvents;
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.Value = ~src.Value;
            end
        end
        
        function findEventsCallback(app, ~, ~)
            % Called when user clicks findEventsBtn. Finds preliminary events based on
            % covariance and peaks/min, plots them, updates them using changepoint
            % algorithm, determines activeEvents based on minDur and minSep, plots final
            % events, plots misc axes, plots ensemble averages, and plots first event on
            % individual event/likelihood axes.
            
            backup = app.storeBackupTab1;
            try
                app.newEventsReset;
                store = app.disable;
                
                % If input had a key, find the true events.
                if app.misc.simInput
                    app.updateCurrentTask('Finding the real events.');
                    key = app.data.origData.Key(app.controls.deselect.indices);
                    app.misc.realEvents = app.findRealEvents(key);
                    app.plotRealEvents(app.misc.realEvents);
                end
                
                % Then find the preliminary events.
                app.updateCurrentTask('Detecting events based on the covariance.');
                cov = app.load.cov.cov;
                switch app.useman
                    case false
                        p1 = app.autop1;
                        p2 = app.autop2;
                        m = app.automin;
                    case true
                        p1 = app.manp1;
                        p2 = app.manp2;
                        m = app.manmin;
                end
                switch app.usemin
                    case false
                        prelimEvents = app.findPrelimEvents(cov, p1, p2);
                    case true
                        prelimEvents = app.findPrelimEventsByMin(cov, m);
                end
                if ~isempty(prelimEvents)
                    % If events exist, plot them.
                    app.plotPrelimEvents(prelimEvents);
                    
                    % Update the events using the changepoint algorithm.
                    finalEvents = app.refineEventsChangepoint(prelimEvents);
                    app.initializeIndividualEvent;
                    app.allEvents = finalEvents;
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end
                    app.data.defaultEvents = app.allEvents;
                    
                    % Allow user to toggle events or save to excel, disable Find Events
                    % button, and show the individual event panel.
                    app.addToNew(true, [app.controls.toggleEventsBtn,...
                        app.excel.addBtn]);
                    app.addToNew(false, app.controls.findEventsBtn);
                    app.analyze.individual.panel.Visible = 'on';
                else
                    % There were no events detected from the covariance.
                    warndlg('No events were found.')
                end
                
                app.controls.findEventsBtn.String = 'Find Events';
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function plotPrelimEvents(app, events)
            % Plots events on the covariance axes in load panel of tab 1.
            
            app.updateCurrentTask('Plotting the initial events.');
            
            % Hide all previous preliminary event patches.
            delete(app.load.cov.patches);
            app.load.cov.patches = gobjects(0);
            
            % Plot new preliminary event patches.
            limy = ylim(app.load.cov.covAxes);
            liminit = limy;
            limy(1) = limy(1) - 0.1*(limy(2) - limy(1)); % Each patch will span from 10% below...
            limy(2) = limy(2) + 0.1*(limy(2) - limy(1)); % ...to 10% above current ylim.
            fileTime = app.data.time;
            for nevent = 1:size(events, 1)
                % For each preliminary event...
                event = events(nevent,:);
                app.load.cov.patches(end+1) = patch(app.load.cov.covAxes,...
                    'XData', [fileTime(event(1)+1) fileTime(event(2)) fileTime(event(2)) fileTime(event(1)+1)],...
                    'YData', [limy(1) limy(1) limy(2) limy(2)],...
                    'FaceColor', app.misc.colors.covPatches,...
                    'FaceAlpha', 0.3,...
                    'HitTest', 'off'); % ...plot a corresponding patch.
            end
            ylim(app.load.cov.covAxes, liminit);
            
            drawnow;
        end
        
        function refinedEvents = refineEventsChangepoint(app, events)
            % Determines start and end times of each event in events using the changepoint
            % algorithm.
            
            app.updateCurrentTask('Determining appropriate window for each event.');
            
            % To store the new events.
            refinedEvents = zeros(size(events));
            
            % Get the beads.
            averageBeads = app.getBeadsToAnalyze;
            
            % Obtain windows for each event.
            [numPtsBefore, numPtsAfter] = app.findChangepointWindows(events, size(averageBeads,1));
            
            % Run the changepoint algorithm on each event.
            if app.controls.useChangepointBtn.Value
                N = size(events, 1);
                for i = 1:N
                    % For each event...
                    app.updateCurrentTask(sprintf('Running the change point algorithm on event %d/%d.', i, N));
                    window = averageBeads(events(i,1) - numPtsBefore(i):events(i,2) + numPtsAfter(i)); % ...consider a window of data surrounding that event...
                    [T1, T2] = app.changepoint(window); % ...and pass this window of data to changepoint().
                    refinedEvents(i,:) = (events(i,1) - numPtsBefore(i) - 1) + [T1, T2]; % Save the resulting changepoints as indices in averageBeads.
                end
            else
                refinedEvents = events;
            end
            
            app.updateCurrentTask('Cleaning up.');
            
            % Store the peaks/min used to determine these events.
            if app.useman
                p1 = app.manp1;
                p2 = app.manp2;
                m = app.manmin;
            else
                p1 = app.autop1;
                p2 = app.autop2;
                m = app.automin;
            end
            if app.usemin
                p1 = [];
                p2 = [];
            else
                m = [];
            end
            app.misc.active.allEvents.cov = app.load.cov.cov;
            app.misc.active.allEvents.peak1 = p1;
            app.misc.active.allEvents.peak2 = p2;
            app.misc.active.allEvents.min = m;
            app.misc.active.allEvents.whichBeads = app.controls.whichBeadBtnGrp.SelectedObject;
            app.misc.active.allEvents.useChangepoint = app.controls.useChangepointBtn.Value;
            
            % The following variables are used for plotting the current event on
            % individual event/likelihood axes.
            app.data.windowStarts = events(:,1) - numPtsBefore;
            app.data.windowStops = events(:,2) + numPtsAfter;
            app.data.defaultWindowStarts = app.data.windowStarts;
            app.data.defaultWindowStops = app.data.windowStops;
        end
        
        function pruneEvents(app)
            % Updates activeEvents based on minDur and minSep. Events shorter than minDur
            % points or separated from other events by fewer than minSep points are
            % 'removed', i.e. the corresponding index in activeEvents is set to false.
            
            app.updateCurrentTask('Removing events by minimum duration and minimum separation.');
            
            dur = app.allEvents(:,2) - app.allEvents(:,1); % Find event durations.
            durcrit = dur >= app.minDur; % Mark which events are at least minDur points long.
            sep = app.allEvents(2:end,1) - app.allEvents(1:end-1,2); % Find event separations.
            sepcrit = logical([1; sep >= app.minSep] & [sep >= app.minSep; 1]); % Mark which events are separated from other events by at least minSep points.
            
            app.analyze.individual.removedByDur = ~durcrit;
            app.analyze.individual.removedBySep = ~sepcrit;
            if isempty(app.analyze.individual.removedByUser)
                % Only update list of which events were removed by user if list is empty,
                % i.e. events were just recently reset. Otherwise, allEvents, minDur, or
                % minSep was just updated, and so this list should not change.
                app.analyze.individual.removedByUser = false(size(app.allEvents,1),1);
            end
            if isempty(app.analyze.individual.includedByUser)
                % Only update list of which events were included by user if list is empty,
                % i.e. events were just recently reset. Otherwise, allEvents, minDur, or
                % minSep was just updated, and so this list should not change.
                app.analyze.individual.includedByUser = false(size(app.allEvents,1),1);
            end
            
            % Mark as valid all events which satisfy both above criteria and were not
            % removed by the user, but mark any event included by user as valid regardless
            % of above criteria.
            app.activeEvents = (durcrit & sepcrit & ~app.analyze.individual.removedByUser)...
                | app.analyze.individual.includedByUser;
            if ~isempty(app.misc.listenerError)
                rethrow(app.misc.listenerError)
            end
        end
        
        function minDurCallback(app, src, ~)
            % Called when user changes the value of minDurEdit. Makes sure input is valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                if ~isfinite(d) || ~isreal(d)
                    % If input is infinite or NaN, reset edit text string to current value and
                    % return.
                    uiwait(errordlg('Minimum duration should be a non-negative integer.', 'Error Dialog', 'modal'))
                    src.String = int2str(app.minDur);
                    app.restore(store);
                    return
                end

                if d < 0
                    % If input is negative, set it to 0.
                    warndlg('Minimum duration should be non-negative. Setting to 0.')
                    d = 0;
                end

                if mod(d,1) ~= 0
                    % If input is a decimal, round down to nearest integer.
                    warndlg('Minimum duration should be an integer. Rounding down.')
                    d = floor(d);
                end

                % Input is a non-negative integer here.
                if app.minDur ~= d
                    % If new value does not match current value, update value and edit text.
                    app.minDur = d;
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end
                    src.String = int2str(d);
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.String = int2str(app.minDur);
            end
        end

        function minSepCallback(app, src, ~)
            % Called when user changes the value of minSepEdit. Makes sure input is valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                if ~isfinite(d) || ~isreal(d)
                    % If input is infinite or NaN, reset edit text string to current value and
                    % return.
                    uiwait(errordlg('Minimum separation should be a non-negative integer.', 'Error Dialog', 'modal'))
                    src.String = int2str(app.minSep);
                    app.restore(store);
                    return
                end

                if d < 0
                    % If input is negative, set it to 0.
                    warndlg('Minimum separation should be non-negative. Setting to 0.')
                    d = 0;
                end

                if mod(d,1) ~= 0
                    % If input is a decimal, round down to nearest integer.
                    warndlg('Minimum separation should be an integer. Rounding down.')
                    d = floor(d);
                end

                % Input is a non-negative integer here.
                if app.minSep ~= d
                    % If new value does not match current value, update value and edit text.
                    app.minSep = d;
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end
                    src.String = int2str(d);
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.String = int2str(app.minSep);
            end
        end
        
        function defcutoffsCallback(app, ~, ~)
            % Called when user clicks cutoffs.resetBtn. Restores minDur and minSep to
            % default values.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Update values.
                app.minDur = app.controls.cutoffs.defMinDur;
                app.minSep = app.controls.cutoffs.defMinSep;
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end

                % Update edit texts.
                app.controls.cutoffs.minDurEdit.String = int2str(app.minDur);
                app.controls.cutoffs.minSepEdit.String = int2str(app.minSep);

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function checkOutdatedEvents(app)
            % Compares the current covariance, peaks, min, and bead(s) to those which were
            % used to determine the present events. If events are out of data, enables a
            % button to let the user recalculate them.
            
            if app.useman
                m = app.manmin;
                p1 = app.manp1;
                p2 = app.manp2;
            else
                m = app.automin;
                p1 = app.autop1;
                p2 = app.autop2;
            end
            if app.usemin
                p1 = [];
                p2 = [];
            else
                m = [];
            end
            
            % Determine if events are outdated...
            covEventsChange = ~isequal(app.load.cov.cov, app.misc.active.allEvents.cov);
            p1EventsChange = ~isequal(p1, app.misc.active.allEvents.peak1);
            p2EventsChange = ~isequal(p2, app.misc.active.allEvents.peak2);
            minEventsChange = ~isequal(m, app.misc.active.allEvents.min);
            beadChange = ~isequal(app.controls.whichBeadBtnGrp.SelectedObject, app.misc.active.allEvents.whichBeads);
            useCPChange = ~isequal(app.controls.useChangepointBtn.Value, app.misc.active.allEvents.useChangepoint);
            enableFindEvents = covEventsChange || minEventsChange || p1EventsChange || p2EventsChange || beadChange || useCPChange;
            
            % ...and update findEventsBtn (and addBtn) accordingly.
            app.addToNew(enableFindEvents, app.controls.findEventsBtn);
            if ~isempty(app.allEvents)
                if enableFindEvents
                    app.controls.findEventsBtn.String = ['<html><p align="center">'...
                        'Recalculate Events<br/>'...
                        'and Update Plots<br/>'...
                        '(this resets change points and<br/>'...
                        'any events removed by user)'... 
                        '</p></html>'];
                    app.addToNew(false, app.excel.addBtn);
                else
                    app.controls.findEventsBtn.String = 'Find Events';
                    app.addToNew(true, app.excel.addBtn);
                end
            end
        end
        
        function plotEvents(app)
            % Plots events on the main axes in load panel of tab 1.
            
            app.updateCurrentTask('Plotting the final events.');
            
            % Hide all previous event patches.
            delete(app.load.main.patches);
            app.load.main.patches = gobjects(0);
            
            % Plot new event patches.
            limy = app.load.main.ylim;
            limy(1) = limy(1) - 0.1*(limy(2) - limy(1)); % Each patch will span from 10% below...
            limy(2) = limy(2) + 0.1*(limy(2) - limy(1)); % ...to 10% above current ylim.
            fileTime = app.data.time;
            for nevent = 1:size(app.allEvents, 1)
                % For each event...
                event = app.allEvents(nevent,:);
                app.load.main.patches(end+1) = patch(app.load.main.mainAxes,...
                    'XData', [fileTime(event(1)+1) fileTime(event(2)) fileTime(event(2)) fileTime(event(1)+1)],...
                    'YData', [limy(1) limy(1) limy(2) limy(2)],...
                    'FaceColor', app.misc.colors.mainPatches,...
                    'FaceAlpha', 0.3,...
                    'Visible', 'off',...
                    'HitTest', 'off'); % ...plot a corresponding patch.
            end
            ylim(app.load.main.mainAxes, app.load.main.ylim);
            
            drawnow;
        end
        
        function updateEventVisNColor(app)
            % Updates Visible and Color properties of event patches to reflect showEvents,
            % activeEvents and current event.
             
            if length(app.load.main.patches) ~= length(app.activeEvents)
                % If the number of event patches do not match the number of events, this
                % function was called before the patches were updated, so hide all patches
                % and return.
                for npatch = 1:length(app.load.main.patches)
                    app.load.main.patches(npatch).Visible = 'off';
                end
                return
            end
            
            j = app.getCurrentEventAsIndex; % Current event.
            for nevent = 1:length(app.activeEvents)
                % For each event...
                if ~app.showEvents
                    % ...if showEvents is false, hide the patch...
                    vis = 'off';
                elseif app.activeEvents(nevent)
                    % ...otherwise, if the event is active, show the patch...
                    vis = 'on';
                elseif nevent == j
                    % ...otherwise, if the event is current, show the patch...
                    vis = 'on';
                else
                    % ...otherwise, event is neither active nor current, so hide the
                    % patch.
                    vis = 'off';
                end
                col = app.misc.colors.mainPatches; % Make the patch red...
                if nevent == j
                    % ...unless it is current, in which case make it gold.
                    col = app.misc.colors.mainPatchesSelected;
                end
                app.load.main.patches(nevent).Visible = vis; % Update Visible.
                app.load.main.patches(nevent).FaceColor = col; % Update Color.
            end
            
            drawnow;
        end
        
        function toggleEventsCallback(app, src, ~)
            % Called when user clicks toggleEventsBtn. Updates showEvents accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.showEvents = src.Value;
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.Value = app.showEvents;
            end
        end
        
        function plotMisc(app)
            % Calculates and plots various data and fits on the misc axes in analyze panel
            % of tab 1. Uses value of app.analyze.misc.data to determine data and, if
            % appropriate, fit. Defaults to event durations if needed.
            
            % Determine which data and distribution to use.
            ylab = 'cumulative distribution';
            switch app.analyze.misc.data
                case 'Step 1 Size'
                    [~, stepA] = app.getPosNStepSizes('A');
                    [~, stepB] = app.getPosNStepSizes('B');
                    if app.misc.simInput
                        [~, stepAr] = app.getRealPosNStepSizes('A');
                        [~, stepBr] = app.getRealPosNStepSizes('B');
                    end
                    if isequal(app.misc.active.allEvents.whichBeads, app.controls.beadABtn)
                        % Just bead A.
                        Y = stepA(:,1);
                        if app.misc.simInput, Yr = stepAr(:,1); end
                        xlab = 'step 1 size (nm), estimated by bead A';
                    elseif isequal(app.misc.active.allEvents.whichBeads, app.controls.beadBBtn)
                        % Just bead B.
                        Y = stepB(:,1);
                        if app.misc.simInput, Yr = stepBr(:,1); end
                        xlab = 'step 1 size (nm), estimated by bead B';
                    else
                        % Both.
                        Y = (stepA(:,1)+stepB(:,1))/2;
                        if app.misc.simInput, Yr = (stepAr(:,1)+stepBr(:,1))/2; end
                        xlab = 'step 1 size (nm), estimated by average of both beads';
                    end
                    type = 'normcdf';
                    app.updateCurrentTask('Plotting the distribution of substep 1 sizes.');
                case 'Step 2 Size'
                    [~, stepA] = app.getPosNStepSizes('A');
                    [~, stepB] = app.getPosNStepSizes('B');
                    if app.misc.simInput
                        [~, stepAr] = app.getRealPosNStepSizes('A');
                        [~, stepBr] = app.getRealPosNStepSizes('B');
                    end
                    if isequal(app.misc.active.allEvents.whichBeads, app.controls.beadABtn)
                        % Just bead A.
                        Y = stepA(:,2);
                        if app.misc.simInput, Yr = stepAr(:,2); end
                        xlab = 'step 2 size (nm), estimated by bead A';
                    elseif isequal(app.misc.active.allEvents.whichBeads, app.controls.beadBBtn)
                        % Just bead B.
                        Y = stepB(:,2);
                        if app.misc.simInput, Yr = stepBr(:,2); end
                        xlab = 'step 2 size (nm), estimated by bead B';
                    else
                        % Both.
                        Y = (stepA(:,2)+stepB(:,2))/2;
                        if app.misc.simInput, Yr = (stepAr(:,2)+stepBr(:,2))/2; end
                        xlab = 'step 2 size (nm), estimated by avg of both beads';
                    end
                    type = 'normcdf';
                    app.updateCurrentTask('Plotting the distribution of substep 2 sizes.');
                case 'Total Step Size'
                    [~, stepA] = app.getPosNStepSizes('A');
                    [~, stepB] = app.getPosNStepSizes('B');
                    if app.misc.simInput
                        [~, stepAr] = app.getRealPosNStepSizes('A');
                        [~, stepBr] = app.getRealPosNStepSizes('B');
                    end
                    if isequal(app.misc.active.allEvents.whichBeads, app.controls.beadABtn)
                        % Just bead A.
                        Y = stepA(:,3);
                        if app.misc.simInput, Yr = stepAr(:,3); end
                        xlab = 'total step size (nm), estimated by bead A';
                    elseif isequal(app.misc.active.allEvents.whichBeads, app.controls.beadBBtn)
                        % Just bead B.
                        Y = stepB(:,3);
                        if app.misc.simInput, Yr = stepBr(:,3); end
                        xlab = 'total step size (nm), estimated by bead B';
                    else
                        % Both.
                        Y = (stepA(:,3)+stepB(:,3))/2;
                        if app.misc.simInput, Yr = (stepAr(:,3)+stepBr(:,3))/2; end
                        xlab = 'total step size (nm), estimated by avg of both beads';
                    end
                    type = 'normcdf';
                    app.updateCurrentTask('Plotting the distribution of total step sizes.');
                case 'Event Duration / Force on Selected Bead'
                    if isequal(app.misc.active.allEvents.whichBeads, app.controls.beadABtn)
                        bead = 'A';
                    else
                        bead = 'B';
                    end
                    force = app.getForce(bead);
                    Y = app.getActiveDur/app.data.fs;
                    X = force(:,2)-force(:,4);
                    type = 'semilogy';
                    ylab = 'event duration (s)';
                    xlab = ['force on bead ' bead ' (pN)'];
                    app.updateCurrentTask(['Plotting the event durations over the force on bead ' bead '.']);
                otherwise
                    Y = app.getActiveDur/app.data.fs;
                    if app.misc.simInput, Yr = app.getRealDur/app.data.fs; end
                    type = 'expcdf';
                    xlab = 'event duration (s)';
                    app.updateCurrentTask('Plotting the distribution of event durations.');
            end
            
            % Determine the data and, if appropriate, the fit to plot.
            switch type
                case 'expcdf'
                    pdf = @(x, k) k*exp(-k*x); % Exponential distribution probability density function.
                    mle = @(k) -sum(log(pdf(Y, k))); % Function to minimize during maximum likelihood estimation.
                    opts = optimset('Display', 'none');
                    [p, ~, exitflag, output] = fminsearch(mle, 0.0005, opts); % k value which minimizes the mle function.
                    if exitflag == 0
                        warndlg(sprintf(['When fitting the event durations to an exponential, fminsearch exited '...
                        'with the following message:\n\n%s\n\nThere may be insufficient data '...
                        'to calculate an accurate fit.'], output.message))
                    end
                    fit = @(x) 1 - exp(-p*x); % Fit to exponential cumulative distribution.
                    n = fit(min(Y))*(length(Y)); % Estimated number of events that are shorter than smallest measured event.
                    X = linspace(min(Y), max(Y), 1000); % Range of D for the cumulative distribution.
                    fitX = linspace(0, max(Y), 1000); % Range of D for the fit.
                    Y = (sum(Y<=X,1)+n) / (length(Y)+n); % Cumulative distribution function.
                    if app.misc.simInput
                        Xr = linspace(min(Yr), max(Yr), 1000);
                        Yr = sum(Yr<=Xr,1) / length(Yr);
                    end
                    fitY = fit(fitX);
                    s = ['k = ' sprintf('%.3g', p) ' s^{-1}'];
                    pos = [0.5 0.5];
                case 'normcdf'
                    pdf = @(x, p) (1/(p(2)*sqrt(2*pi)))*exp(-0.5*((x-p(1))/p(2)).^2); % Normal distribution probability density function.
                    mle = @(p) -sum(log(pdf(Y, p))); % Function to minimize during maximum likelihood estimation.
                    opts = optimset('Display', 'none');
                    [p, ~, exitflag, output] = fminsearch(mle, [mean(Y) std(Y)], opts); % p value which minimizes the mle function.
                    if exitflag == 0
                        warndlg(sprintf(['When fitting the step sizes to a Gaussian, fminsearch exited '...
                            'with the following message:\n\n%s\n\nThere may be insufficient data '...
                            'to calculate an accurate fit.'], output.message))
                    end
                    fit = @(x) 0.5*erfc((p(1)-x)/(p(2)*sqrt(2))); % Fit to normal cumulative distribution.
                    X = linspace(min(Y), max(Y), 1000); % Range of D for the cumulative distribution.
                    fitX = p(1)+p(2)*-sqrt(2)*erfcinv(2*linspace(0.001,0.999,1000)); % Range of D for the fit.
                    Y = sum(Y<=X,1) / length(Y); % Cumulative distribution function.
                    if app.misc.simInput
                        Xr = linspace(min(Yr), max(Yr), 1000);
                        Yr = sum(Yr<=Xr,1) / length(Yr);
                    end
                    fitY = fit(fitX);
                    s = sprintf('\\mu = %.3g nm\\newline\\sigma = %.3g nm', [p(1) p(2)]);
                    pos = [0.25, 5];
                case 'semilogy'
                    fitX = [];
                    fitY = [];
                    s = '';
                    pos = [0, 0];
            end
            
            if app.misc.simInput && (strcmp(type, 'normcdf') || strcmp(type, 'expcdf'))
                % If user is analyzing simulated data and looking at the event durations or step
                % sizes, show the actual distribution in yellow.
                app.changeData(app.analyze.misc.realPlot, Xr, Yr);
            else
                app.changeData(app.analyze.misc.realPlot, 0, 0);
            end
            
            % Update plots.
            app.changeData(app.analyze.misc.miscPlot, X, Y);
            app.changeData(app.analyze.misc.miscFitPlot, fitX, fitY);
            xlabel(app.analyze.misc.miscAxes, xlab);
            ylabel(app.analyze.misc.miscAxes, ylab);
            
            % If scatter plot, change marker, line style, y scale, axis limits.
            if strcmp(type, 'semilogy')
                app.analyze.misc.miscPlot.Marker = 'o';
                app.analyze.misc.miscPlot.LineStyle = 'none';
                app.analyze.misc.miscAxes.YScale = 'log';
                axis(app.analyze2.misc.miscAxes, 'tight');
                limx = app.analyze.misc.miscAxes.XLim;
                limy = app.analyze.misc.miscAxes.YLim;
                app.analyze.misc.miscAxes.XLim = limx + [-0.1 0.1]*(limx(2)-limx(1));
                app.analyze.misc.miscAxes.YLim = 10.^(log10(limy) + [-0.1 0.1]*(log10(limy(2))-log10(limy(1))));
            else
                app.analyze.misc.miscPlot.Marker = 'none';
                app.analyze.misc.miscPlot.LineStyle = '-';
                app.analyze.misc.miscAxes.YScale = 'linear';
                axis(app.analyze.misc.miscAxes, 'tight');
            end
            app.resetAxesDefaultLimits(app.analyze.misc.miscAxes);
            
            % Update label.
            app.analyze.misc.miscLabel.String = '';
            [~, x, y] = app.placeText(app.analyze.misc.miscAxes, s, pos, app.misc.fontSize);
            app.analyze.misc.miscLabel.Position = [x y];
            app.analyze.misc.miscLabel.String = s;
            app.analyze.misc.miscLabel.Visible = 'on';
            
            % Show axes and toolbar.
            app.analyze.misc.miscAxes.Visible = 'on';
            
            % Resize axes based on yticklabels and ylabel.
            ti = app.analyze.misc.miscAxes.TightInset;
            pos = app.analyze.misc.miscAxes.Position;
            pos(3) = 0.66 - pos(1) - ti(3);
            app.analyze.misc.miscAxes.Position = pos;
            
            % Update stored parameters.
            switch app.analyze.misc.data
                case 'Step 1 Size'
                    app.analyze.misc.meanStep1 = p(1);
                    app.analyze.misc.stdStep1 = p(2);
                case 'Step 2 Size'
                    app.analyze.misc.meanStep2 = p(1);
                    app.analyze.misc.stdStep2 = p(2);
                case 'Total Step Size'
                    app.analyze.misc.meanTotalStep = p(1);
                    app.analyze.misc.stdTotalStep = p(2);
                case 'Event Duration / Force on Selected Bead'
                otherwise
                    app.analyze.misc.kDur = p;
            end
            
            drawnow;
        end
        
        function [avgbeg, avgend] = calcEnsemble(app)
            % Calculates the ensemble averages.
            
            app.updateCurrentTask('Calculating the ensemble averages.');
            
            % Get beads, events, and durations.
            averageBeads = app.getBeadsToAnalyze;
            events = app.getActiveEvents;
            dur = app.getActiveDur;
            maxDur = max(dur);
            
            pts_bef =  round(app.analyze.ensemble.pts_bef * app.data.fs);
            pts_after = round(app.analyze.ensemble.pts_after * app.data.fs);
            exten_pos = round(app.analyze.ensemble.exten_pos * app.data.fs);
            
            % Calculate average of event beginnings.
            app.updateCurrentTask('Calculating the event start ensemble average.');
            beginning = ones(maxDur + pts_bef - exten_pos, size(events, 1));
            for i = 1:size(events, 1)
                % For each event...
                if events(i,1) - pts_bef < 1
                    % ...if event window starts N points before start of file, then first
                    % N values of event window will match first value in file...
                    event = [averageBeads(1)*ones(pts_bef - events(i,1),1); ...
                        averageBeads(1:events(i,2) - exten_pos)];
                else
                    % ...otherwise determine the event window...
                    event = averageBeads(events(i,1) + 1 - pts_bef:events(i,2) - exten_pos);
                end
                % ...and calculate the extension so that each event beginning is the same
                % size as other event beginnings.
                extension = mean(averageBeads(events(i,2) - exten_pos - app.data.fs/200:events(i,2) - exten_pos))*ones(maxDur - dur(i),1);
                beginning(:,i) = [event; extension];
            end
            avgbeg = mean(beginning, 2); % Average the event beginnings.
            
            % Calculate average of event endings.
            app.updateCurrentTask('Calculating the event end ensemble average.');
            ending = ones(maxDur + pts_after - exten_pos, size(events, 1));
            for i = 1:size(events, 1)
                % For each event...
                if events(i,2) + pts_after > length(averageBeads)
                    % ...if event window ends N points after end of file, then last N
                    % values of event window will match last value in file...
                    event = [averageBeads(events(i,1) + 1 + exten_pos:end); ...
                        averageBeads(end)*ones(pts_after + events(i,2) - length(averageBeads),1)];
                else
                    % ...otherwise determine the event window...
                    event = averageBeads(events(i,1) + 1 + exten_pos:events(i,2) + pts_after);
                end
                % ...and calculate the extension so that each event ending is the same
                % size as other event endings.
                extension = mean(averageBeads(events(i,1) + 1 + exten_pos:events(i,1) + 1 + exten_pos + app.data.fs/200))*ones(maxDur - dur(i),1);
                ending(:,i) = [extension; event];
            end
            avgend = mean(ending, 2); % Average the event endings.
            
            % Adjust averages so beginning average starts at 0.
            app.updateCurrentTask('Cleaning up.');
            offset = movmean(avgbeg(1:pts_bef), app.data.fs/10);
            avgbeg = avgbeg - offset(end);
            avgend = avgend - offset(end);
        end
        
        function plotEnsemble(app)
            % Plots the ensemble averages on the ensemble axes in analyze panel of tab 1.
            
            app.updateCurrentTask('Plotting the ensemble averages.');
            
            startAvg = app.analyze.ensemble.start;
            endAvg = app.analyze.ensemble.end;
            if app.flipped
                % If data should be flipped, then averages should be flipped.
                startAvg = -startAvg;
                endAvg = -endAvg;
            end
            
            if isequal(app.misc.active.allEvents.whichBeads, app.controls.beadABtn)
            	ylab = 'position of bead A (nm)';
            elseif isequal(app.misc.active.allEvents.whichBeads, app.controls.beadBBtn)
                ylab = 'position of bead B (nm)';
            else
                ylab = 'avg pos of both beads (nm)';
            end
            
            % Calculate time ranges for each average.
            startTime = (1:length(startAvg)) / app.data.fs;
            endTime = (length(startTime)+(1:length(endAvg))) / app.data.fs;
            
            % Update plots.
            app.changeData(app.analyze.ensemble.startPlot, startTime, startAvg);
            app.changeData(app.analyze.ensemble.endPlot, endTime, endAvg);
            axis(app.analyze.ensemble.ensembleAxes, 'tight');
            app.analyze.ensemble.ensembleAxes.Visible = 'on';
            ylabel(app.analyze.ensemble.ensembleAxes, ylab);
            
            % Resize axes based on yticklabels and ylabel.
            ti = app.analyze.ensemble.ensembleAxes.TightInset;
            pos = app.analyze.ensemble.ensembleAxes.Position;
            pos(3) = 0.99 - pos(1) - ti(3);
            app.analyze.ensemble.ensembleAxes.Position = pos;
            
            % Hide fits and labels.
            app.analyze.ensemble.startFitPlot.Visible = 'off';
            app.analyze.ensemble.endFitPlot.Visible = 'off';
            app.analyze.ensemble.startLabel.Visible = 'off';
            app.analyze.ensemble.endLabel.Visible = 'off';
            app.analyze.ensemble.stepLabel.Visible = 'off';
            
            drawnow;
        end
        
        function [startFit, endFit, p] = calcEnsembleFits(app)
            % Calculates fits of the ensemble averages.
            
            app.updateCurrentTask('Fitting the ensemble averages.');
            
            fs = app.data.fs;
            pts_bef = round(app.analyze.ensemble.pts_bef * fs);
            pts_after = round(app.analyze.ensemble.pts_after * fs);
            skip = round(app.analyze.ensemble.num_pts_to_skip);
            
            % Isolate the exponential portions of each average.
            avgbegexp = app.analyze.ensemble.start(pts_bef+1+skip:end);
            startTime = (skip-1+(1:length(avgbegexp))) / fs;
            avgendexp = app.analyze.ensemble.end(1:end-pts_after-skip);
            endTime = ((-length(avgendexp):-1)+1-skip) / fs;
            
            % Initialize variables.
            p = [];
            startFit = NaN(1, length(avgbegexp));
            endFit = NaN(1, length(avgendexp));
            t = timer('StartFcn', @(varargin) set(app.load.skipBtn, 'Enable', 'on'),...
                'TimerFcn', @(varargin) set(app.load.skipBtn, 'Visible', 'on'),...
                'StartDelay', 3);
            
            app.updateCurrentTask('Fitting the ensemble averages.');
            try
                start(t)
                [startFit, endFit, p] = app.fitexp(avgbegexp, startTime, avgendexp, endTime, app.analyze.ensemble.globalFit); % [start fit, end fit, fit parameters]
                stop(t)
                app.load.skipBtn.Visible = 'off';
            catch
                stop(t)
                app.load.skipBtn.Visible = 'off';
            end
        end
        
        function [fstart, fend, p] = fitexp(app, ystart, xstart, yend, xend, globalFit)
            % Fits exponentials to data.
            
            % Inputs should be column vectors.
            ystart = ystart(:);
            yend = yend(:);
            xstart = xstart(:);
            xend = xend(:);

            nstart = length(ystart);
            
            % Fit using lsqcurvefit.
            % REQUIRES OPTIMIZATION TOOLBOX
            options = optimoptions('lsqcurvefit',...
                'MaxFunctionEvaluations', 100000,...
                'MaxIterations', 100000,...
                'FunctionTolerance', 10e-12,...
                'Display', 'off');
            app.misc.skip = false; % To be sure that a previous skip button click won't affect this upcoming fit.
            
            if globalFit
                fun = @(p,x) expFunGlobal(p,x,nstart);
                
                % Determine initial guesses and parameter bounds, which depend on expected
                % shape of data.
                p0 = [yend(end)-yend(1), 2, 2, yend(1)]; % [a kstart kend y0]
                if mean(yend) > 0 % Upright.
                    lb = [0, 0, 0, 0];
                    ub = [Inf, Inf, Inf, Inf];
                else % Upside down.
                    lb = [-Inf, 0, 0, -Inf];
                    ub = [0, Inf, Inf, 0];
                end
                
                t = tic;
                p = lsqcurvefit(fun, p0, [xstart; xend], [ystart; yend], lb, ub, options);
                kstart = p(2);
                kend = p(3);
                step1 = p(4);
                totalstep = p(1)+p(4);

                % Prepare output.
                f = fun(p, [xstart; xend]);
                fstart = f(1:nstart);
                fend = f(nstart+1:end);
                
                p = [kstart kend step1 totalstep];
            else
                fun = @expFun;

                % Determine initial guesses and parameter bounds, which depend on expected
                % shape of data.
                p0 = [ystart(1)-ystart(end), -2, ystart(end)]; % [a kstart y0]
                if mean(ystart) > 0 % Upright.
                    lb = [-Inf -Inf 0];
                    ub = [0 0 Inf];
                else % Upside down.
                    lb = [0 -Inf -Inf];
                    ub = [Inf 0 0];
                end
                
                t = tic;
                p = lsqcurvefit(fun, p0, xstart, ystart, lb, ub, options);
                fstart = fun(p, xstart);
                kstart = -p(2);
                totalstep = p(3);
                
                p0 = [yend(end)-yend(1), 2, yend(1)]; % [a kend y0]
                if mean(yend) > 0 % Upright.
                    lb = [0 0 0];
                    ub = [Inf Inf Inf];
                else % Upside down.
                    lb = [-Inf 0 -Inf];
                    ub = [0 Inf 0];
                end
                
                t = tic;
                p = lsqcurvefit(fun, p0, xend, yend, lb, ub, options);
                fend = fun(p, xend);
                kend = p(2);
                step1 = p(3);
                
                p = [kstart kend step1 totalstep];
            end
            
            function val = expFunGlobal(p, x, n)
                persistent count;
                if isempty(count), count = 0; end
                t0 = toc(t);
                c0 = floor(t0/1); % Call drawnow every second.
                if ~isequal(count, c0)
                    count = c0;
                    drawnow;
                    if app.misc.skip, error('expFun:abort', 'User clicked skip button.'); end
                end
                val = [-p(1)*exp(-p(2)*x(1:n)) + (p(4)+p(1));
                    p(1)*exp(p(3)*x(n+1:end)) + p(4)];
            end
            
            function val = expFun(p, x)
                persistent count;
                if isempty(count), count = 0; end
                t0 = toc(t);
                c0 = floor(t0/1); % Call drawnow every second.
                if ~isequal(count, c0)
                    count = c0;
                    drawnow;
                    if app.misc.skip, error('expFun:abort', 'User clicked skip button.'); end
                end
                val = p(1)*exp(p(2)*x) + p(3);
            end
        end
        
        function plotEnsembleFits(app)
            % Plots the ensemble average fits and labels on the ensemble axes in analyze
            % panel of tab 1.
            
            app.updateCurrentTask('Plotting the ensemble average fits.');
            
            startFit = app.analyze.ensemble.startFit;
            endFit = app.analyze.ensemble.endFit;
            if app.flipped
                % If data should be flipped, then fits should be flipped.
                startFit = -startFit;
                endFit = -endFit;
            end
            pts_bef = round(app.analyze.ensemble.pts_bef * app.data.fs);
            skip = round(app.analyze.ensemble.num_pts_to_skip);
            
            % Erase all labels so they don't mess up new label positions (they should
            % already be invisible due to plotEnsemble()).
            app.analyze.ensemble.startLabel.String = '';
            app.analyze.ensemble.endLabel.String = '';
            app.analyze.ensemble.stepLabel.String = '';
            
            if ~any(isnan(startFit)) && ~isempty(app.analyze.ensemble.p)
                % If start fit exists, update the start fit plot...
                startTime = (pts_bef+1+skip+(1:length(startFit))) / app.data.fs;
                app.changeData(app.analyze.ensemble.startFitPlot, startTime, startFit);
                
                % ...and label.
                kf = app.analyze.ensemble.p(1);
                s = sprintf('k_f = %.3g s^{-1}', kf);
                [~, x, y] = app.placeText(app.analyze.ensemble.ensembleAxes, s, [0.25, 0.5], app.misc.fontSize);
                app.analyze.ensemble.startLabel.Position = [x y];
                app.analyze.ensemble.startLabel.String = s;
                app.analyze.ensemble.startLabel.Visible = 'on';
            end
            
            if ~any(isnan(endFit)) && ~isempty(app.analyze.ensemble.p)
                % If end fit exists, update the end fit plot...
                endTime = (pts_bef+1+skip+length(startFit)+(1:length(endFit))) / app.data.fs;
                app.changeData(app.analyze.ensemble.endFitPlot, endTime, endFit);
                
                % ...and label.
                kr = app.analyze.ensemble.p(2);
                s = sprintf('k_r = %.3g s^{-1}', kr);
                [~, x, y] = app.placeText(app.analyze.ensemble.ensembleAxes, s, [0.75, 0.5], app.misc.fontSize);
                app.analyze.ensemble.endLabel.Position = [x y];
                app.analyze.ensemble.endLabel.String = s;
                app.analyze.ensemble.endLabel.Visible = 'on';
            end
            
            if ~isempty(app.analyze.ensemble.p)
                % Update step labels.
                step1 = app.analyze.ensemble.p(3);
                totalStep = app.analyze.ensemble.p(4);
                if app.flipped
                    step1 = -step1;
                    totalStep = -totalStep;
                end
                s = sprintf('step 1: %.3g nm\ntotal step: %.3g nm', step1, totalStep);
                limy = ylim(app.analyze.ensemble.ensembleAxes);
                [~, x, y] = app.placeText(app.analyze.ensemble.ensembleAxes, s, [0.5, -limy(1)/(limy(2)-limy(1))], app.misc.fontSize);
                app.analyze.ensemble.stepLabel.Position = [x y];
                app.analyze.ensemble.stepLabel.String = s;
                app.analyze.ensemble.stepLabel.Visible = 'on';
            end
            
            axis(app.analyze.ensemble.ensembleAxes, 'tight');
            
            drawnow;
        end
        
        function initializeIndividualEvent(app)
            % Plots all data on individual event axes. Plots true events if data is
            % simulated.
            
            averageBeads = app.getBeadsToAnalyze;
            if app.flipped
                % If data should be flipped, flip it.
                averageBeads = -averageBeads;
            end
            fileTime = app.data.time;
            app.changeData(app.analyze.individual.eventPlot, fileTime, averageBeads);
            app.analyze.individual.eventPlot.Visible = 'off';
            
            % If file is simulated, plot all true events within the window.
            if app.misc.simInput
                app.updateCurrentTask('Plotting real events');
                
                % Initialize lines.
                delete(app.misc.simLines);
                app.misc.simLines = gobjects(0);
                
                % For each true event start/end, plot a yellow vertical line.
                limy = [min(averageBeads) max(averageBeads)];
                for x = app.misc.realEvents(:).'
                    app.misc.simLines(end+1) = plot(app.analyze.individual.eventAxes,...
                        fileTime([x x]), limy,...
                        'Color', app.misc.colors.simLines,...
                        'HitTest', 'off',...
                        'Visible', 'off');
                end
            end
            
            if isequal(app.misc.active.allEvents.whichBeads, app.controls.beadABtn)
            	ylab = 'position of bead A (nm)';
            elseif isequal(app.misc.active.allEvents.whichBeads, app.controls.beadBBtn)
                ylab = 'position of bead B (nm)';
            else
                ylab = 'avg pos of both beads (nm)';
            end
            ylabel(app.analyze.individual.eventAxes, ylab);
        end
        
        function plotIndividualEvent(app)
            % Plots the current event on the individual event/likelihood axes in analyze
            % panel of tab 1.
            
            delete(getappdata(app.analyze.individual.eventAxes, 'listener'))
            
            % Get events, current event, beads, and time.
            events = app.allEvents;
            j = app.getCurrentEventAsIndex;
            if isempty(j) || j < 1 || j > size(events, 1)
                % Exit if current event does not exist or is not within valid range.
                return
            end
            averageBeads = app.analyze.individual.eventPlot.YData;
            fileTime = app.analyze.individual.eventPlot.XData;
            
            % Clear individual event axes so y limit is based only on plot of event j.
            app.changeData(app.analyze.individual.defTonPlot, [], []);
            app.changeData(app.analyze.individual.defToffPlot, [], []);
            app.changeData(app.analyze.individual.defMean1Plot, [], []);
            app.changeData(app.analyze.individual.defUpLim1Plot, [], []);
            app.changeData(app.analyze.individual.defLowLim1Plot, [], []);
            app.changeData(app.analyze.individual.defMean2Plot, [], []);
            app.changeData(app.analyze.individual.defUpLim2Plot, [], []);
            app.changeData(app.analyze.individual.defLowLim2Plot, [], []);
            app.changeData(app.analyze.individual.defMean3Plot, [], []);
            app.changeData(app.analyze.individual.defUpLim3Plot, [], []);
            app.changeData(app.analyze.individual.defLowLim3Plot, [], []);
            
            app.updateCurrentTask(sprintf('Preparing to plot event #%d.', app.analyze.individual.currentEvent));
            
            % Get current event's window, changepoints, and likelihood.
            start = app.data.windowStarts(j);
            stop = app.data.windowStops(j);
            eventTime = fileTime(start:stop); % Individual event plot X data.
            eventData = averageBeads(start:stop); % Individual event plot Y data.
            N = numel(eventData); % Number of points in event window.
            eventTon = events(j,1)-start+1; % First changepoint as index in eventTime.
            eventToff = events(j,2)-start+1; % Second changepoint as index in eventTime.
            [~, ~, Lr, idx] = app.changepoint(eventData, 'likelihood'); % [likelihood, indices of eventTime for which likelihood has been calculated]
            M = size(Lr, 1); % Number of points considered in likelihood calculation.
            [~, LrTon] = min(abs(idx-eventTon)); % First changepoint as index in Lr.
            [~, LrToff] = min(abs(idx-eventToff)); % Second changepoint as index in Lr.
            Z = Lr+Lr'; % Likelihood plot Z data.
            
            % Calculate the means and standard deviations of the event window (1) before
            % the first changepoint, (2) between the two changepoints, and (3) after the
            % second changepoint.
            y = cumsum(eventData);
            y2 = cumsum(eventData.^2);
            x2 = (y(eventToff)-y(eventTon))/(eventToff-eventTon);
            x13 = (y(N)-y(eventToff)+y(eventTon))/(N-eventToff+eventTon);
            s213 = ((y2(N)-y2(eventToff)+y2(eventTon)) - ((y(N)-y(eventToff)+y(eventTon)).^2./(N-eventToff+eventTon)))./(N-eventToff+eventTon);
            s213(s213 < 0) = 0;
            s22 = ((y2(eventToff)-y2(eventTon))-((y(eventToff)-y(eventTon))^2/(eventToff-eventTon)))/(eventToff-eventTon);
            s22(s22 < 0) = 0;
            
            % Zoom in on event j and determine y limit.
            app.updateCurrentTask(sprintf('Plotting event #%d.', app.analyze.individual.currentEvent));
            app.analyze.individual.eventPlot.Visible = 'on';
            if app.misc.simInput
                for nline = 1:numel(app.misc.simLines)
                    if strcmp(app.misc.simLines(nline).Visible, 'on'), break, end
                    app.misc.simLines(nline).Visible = 'on';
                end
            end
            limy = [min(averageBeads), max(averageBeads)];
            app.analyze.individual.eventAxes.XLim = fileTime([start stop]);
            app.resetAxesDefaultLimits(app.analyze.individual.eventAxes);
            
            % Plot the changepoints, means, and variances on top of beads.
            app.updateCurrentTask(sprintf('Plotting the change points on top of event #%d.', app.analyze.individual.currentEvent));
            app.changeData(app.analyze.individual.defTonPlot, [eventTime(eventTon) eventTime(eventTon)], limy);
            app.changeData(app.analyze.individual.defToffPlot, [eventTime(eventToff) eventTime(eventToff)], limy);
            uistack(app.analyze.individual.defTonPlot, 'top');
            uistack(app.analyze.individual.defToffPlot, 'top');
            app.changeData(app.analyze.individual.defMean1Plot, eventTime(1:eventTon), x13*ones(eventTon,1));
            app.changeData(app.analyze.individual.defUpLim1Plot, eventTime(1:eventTon), (x13+sqrt(s213))*ones(eventTon,1));
            app.changeData(app.analyze.individual.defLowLim1Plot, eventTime(1:eventTon), (x13-sqrt(s213))*ones(eventTon,1));
            app.changeData(app.analyze.individual.defMean2Plot, eventTime(eventTon:eventToff), x2*ones(eventToff-eventTon+1,1));
            app.changeData(app.analyze.individual.defUpLim2Plot, eventTime(eventTon:eventToff), (x2+sqrt(s22))*ones(eventToff-eventTon+1,1));
            app.changeData(app.analyze.individual.defLowLim2Plot, eventTime(eventTon:eventToff), (x2-sqrt(s22))*ones(eventToff-eventTon+1,1));
            app.changeData(app.analyze.individual.defMean3Plot, eventTime(eventToff:end), x13*ones(N-eventToff+1,1));
            app.changeData(app.analyze.individual.defUpLim3Plot, eventTime(eventToff:end), (x13+sqrt(s213))*ones(N-eventToff+1,1));
            app.changeData(app.analyze.individual.defLowLim3Plot, eventTime(eventToff:end), (x13-sqrt(s213))*ones(N-eventToff+1,1));
            
            % Finalize individual event axes.
            ylim(app.analyze.individual.eventAxes, limy);
            app.analyze.individual.eventAxes.Visible = 'on';
            
            % Update plots on likelihood axes.
            app.updateCurrentTask(sprintf('Plotting the likelihood for event #%d.', app.analyze.individual.currentEvent));
            app.changeData(app.analyze.individual.surfPlot, eventTime(idx), eventTime(idx), Z);
            view(app.analyze.individual.likelihoodAxes, 2);
            app.changeData(app.analyze.individual.defXonPlot, eventTime(idx), eventTime(eventTon)*ones(M,1), Z(LrTon,:));
            app.changeData(app.analyze.individual.defXoffPlot, eventTime(idx), eventTime(eventToff)*ones(M,1), Z(LrToff,:));
            app.changeData(app.analyze.individual.defYonPlot, eventTime(eventTon)*ones(M,1), eventTime(idx), Z(LrTon,:));
            app.changeData(app.analyze.individual.defYoffPlot, eventTime(eventToff)*ones(M,1), eventTime(idx), Z(LrToff,:));
            
            % Finalize likelihood axes.
            axis(app.analyze.individual.likelihoodAxes, 'tight');
            app.resetAxesDefaultLimits(app.analyze.individual.likelihoodAxes);
            app.analyze.individual.likelihoodAxes.Visible = 'on';
            lh = addlistener(app.analyze.individual.eventAxes, 'XLim', 'PostSet', @app.eventXAxisUpdate);
            setappdata(app.analyze.individual.eventAxes, 'listener', lh);
            
            % Resize axes based on yticklabels and ylabel.
            app.updateCurrentTask('Cleaning up.');
            tiL = app.analyze.individual.likelihoodAxes.TightInset;
            tiE = app.analyze.individual.eventAxes.TightInset;
            posL = app.analyze.individual.likelihoodAxes.Position;
            posE = app.analyze.individual.eventAxes.Position;
            posE(1) = posL(1) + posL(3) + tiL(3) + 0.01;
            posE(3) = 0.745 - posE(1) - tiE(3);
            app.analyze.individual.eventAxes.Position = posE;
            
            % Show buttons that allow interaction with likelihood plot.
            app.analyze.individual.pickChngptBtn.Visible = 'on';
            app.analyze.individual.resetBtn.Visible = 'on';
            app.analyze.individual.snap.Visible = 'on';
            
            drawnow;
        end
        
        function eventXAxisUpdate(app, ~, ~)
            % Called when user pans individual event plot. Calculates the likelihood for
            % whatever new window of data is showing.
            
            try
                % Calculate new likelihood and nearest likelihood values for the current
                % changepoints.
                fileTime = app.analyze.individual.eventPlot.XData;
                averageBeads = app.analyze.individual.eventPlot.YData;
                fileStart = app.analyze.individual.eventAxes.XLim(1); % Time of window start.
                fileStop = app.analyze.individual.eventAxes.XLim(2); % Time of window end.
                [~, start] = min(abs(fileTime-fileStart)); % Window start as index in fileTime.
                [~, stop] = min(abs(fileTime-fileStop)); % Window end as index in fileTime.
                eventTime = fileTime(start:stop);
                eventData = averageBeads(start:stop);
                [~, ~, Lr, idx] = app.changepoint(eventData, 'likelihood'); % [likelihood, indices of eventTime for which likelihood has been calculated]
                Z = Lr+Lr';
                M = size(Lr, 1); % Number of points considered in likelihood calculation.
                events = app.allEvents;
                j = app.getCurrentEventAsIndex;
                eventTon = events(j,1)-start+1; % First changepoint as index in eventTime.
                eventToff = events(j,2)-start+1; % Second changepoint as index in eventTime.
                [~, LrTon] = min(abs(idx-eventTon)); % First changepoint as index in Lr.
                [~, LrToff] = min(abs(idx-eventToff)); % Second changepoint as index in Lr.
                
                % Update likelihood plot.
                app.analyze.individual.likelihoodAxes.XLim = app.analyze.individual.eventAxes.XLim;
                app.analyze.individual.likelihoodAxes.YLim = app.analyze.individual.eventAxes.XLim;
                app.resetAxesDefaultLimits(app.analyze.individual.likelihoodAxes);
                app.changeData(app.analyze.individual.surfPlot, eventTime(idx), eventTime(idx), Z);
                if eventTon > 0 && eventTon < numel(eventTime)
                    app.changeData(app.analyze.individual.defXonPlot, eventTime(idx), eventTime(eventTon)*ones(M,1), Z(LrTon,:));
                    app.changeData(app.analyze.individual.defYonPlot, eventTime(eventTon)*ones(M,1), eventTime(idx), Z(LrTon,:));
                else
                    app.changeData(app.analyze.individual.defXonPlot, [], [], []);
                    app.changeData(app.analyze.individual.defYonPlot, [], [], []);
                end
                if eventToff > 0 && eventToff < numel(eventTime)
                    app.changeData(app.analyze.individual.defXoffPlot, eventTime(idx), eventTime(eventToff)*ones(M,1), Z(LrToff,:));
                    app.changeData(app.analyze.individual.defYoffPlot, eventTime(eventToff)*ones(M,1), eventTime(idx), Z(LrToff,:));
                else
                    app.changeData(app.analyze.individual.defXoffPlot, [], [], []);
                    app.changeData(app.analyze.individual.defYoffPlot, [], [], []);
                end
            catch
                if ~isvalid(app)
                    return
                end
                % Handling errors in this callback is tricky because the user might be
                % panning, e.g., and turning pan off while the mouse is clicked seems to
                % mess things up. We should not interrupt MATLAB from appropriately
                % reacting to further user interaction, such as unclicking the mouse. So,
                % the likelihood plot may be out of date if error occurs, but (1) an error
                % shouldn't occur, (2) if an error does occur, the plot should refresh
                % once this function is called again, provided it doesn't error again, and
                % (3) at the very least, plotIndividualEvent() will update the plot after
                % user changes the current event.
            end
        end
        
        function pickChngptCallback(app, ~, ~)
            % Called when user clicks pickChngptBtn. Allows user to then click on the
            % likelihood axes or, alternatively, click twice on the individual event axes
            % to manually set the changepoints.
            
            backup = app.storeBackupTab1;
            try
                % Temporarily disable click callback so user doesn't click before store is
                % passed to the callback. It will be set once the user has moused over
                % either plot.
                app.misc.storeClickCallback = app.misc.fig.WindowButtonDownFcn;
                app.misc.fig.WindowButtonDownFcn = [];

                store = app.disable;
                drawnow;
                app.analyze.individual.first = true;
                app.misc.storeMotionCallback = app.misc.fig.WindowButtonMotionFcn;
                app.misc.fig.WindowButtonMotionFcn = {@app.pickChngptMotion store backup}; % Set the figure motion callback to pickChngptMotion().
                app.misc.fig.WindowKeyPressFcn = {@app.pickChngptEscape store backup}; % Let user escape by pressing escape key.
                app.updateCurrentTask('User is setting the change points.'); % Update the current task message;

                app.analyze.individual.likelihoodAxes.XColor = app.misc.colors.activeAxes;
                app.analyze.individual.likelihoodAxes.YColor = app.misc.colors.activeAxes;
                app.analyze.individual.likelihoodAxes.ZColor = app.misc.colors.activeAxes;
                app.analyze.individual.eventAxes.XColor = app.misc.colors.activeAxes;
                app.analyze.individual.eventAxes.YColor = app.misc.colors.activeAxes;
                app.analyze.individual.eventAxes.ZColor = app.misc.colors.activeAxes;
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function defChngptCallback(app, ~, ~)
            % Called when user clicks individual.resetBtn. Sets the changepoints for the
            % current event to their original values and updates plots accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                
                % Update allEvents to reflect new changepoints.
                j = app.getCurrentEventAsIndex;
                app.data.windowStarts(j) = app.data.defaultWindowStarts(j);
                app.data.windowStops(j) = app.data.defaultWindowStops(j);
                app.allEvents(j,:) = app.data.defaultEvents(j,:);
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                
                % Update plot.
                app.plotIndividualEvent;
                
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
                
        function pickChngptMotion(app, ~, ~, store, backup)
            % Motion callback for the figure when user is setting manual changepoints by
            % clicking the likelihood plot or individual event plot.
            
            thisMotionCallback = app.misc.fig.WindowButtonMotionFcn;
            try
                if ~isequal(app.misc.tab1, app.misc.tabgroup.SelectedTab)
                    return
                end
                
                fileTime = app.analyze.individual.eventPlot.XData;
                averageBeads = app.analyze.individual.eventPlot.YData;
                fileStart = app.analyze.individual.eventAxes.XLim(1);
                fileStop = app.analyze.individual.eventAxes.XLim(2);
                [~, start] = min(abs(fileTime-fileStart));
                [~, stop] = min(abs(fileTime-fileStop));
                eventTime = fileTime(start:stop);
                eventData = averageBeads(start:stop);
                N = numel(eventData); % Number of points in event window.
                X = app.analyze.individual.surfPlot.XData;
                Y = app.analyze.individual.surfPlot.YData;
                Z = app.analyze.individual.surfPlot.ZData;
                M = size(Z, 1);
                
                [in, p1, p2] = app.inAxes(...                 % [whether in either axes,
                    [app.analyze.individual.likelihoodAxes... % front mouse position in axes,
                    app.analyze.individual.eventAxes]);       % back mouse position in axes]
                if in(1)
                    % If mouse is in the likelihood axes:
                    
                    % Get x and y values (i.e. times) of point nearest to cursor.
                    [x, y] = app.findClosestPoint(p1(1,:), p2(1,:), X, Y, Z);
                    
                    if ~isempty(x)
                        % If user is hovering over the likelihood surface plot:
                        
                        % Determine the new changepoints based on the mouse position.
                        Ton = min([x y]); % Time of first changepoint.
                        Toff = max([x y]); % Time of second changepoint.
                        if app.analyze.individual.snap.Value
                            % If user wants to use snap to changepoint feature:
                            [~, eventTon] = min(abs(eventTime-Ton)); % User-defined first changepoint as index in eventTime.
                            [~, eventToff] = min(abs(eventTime-Toff)); % User-defined second changepoint as index in eventTime.
                            onWindowStart = max(eventTon-25,1); % Start of window for Ton, over which likelihood will be calculated.
                            onWindowData = eventData(onWindowStart:min(eventTon+25,N));
                            onWindowTime = eventTime(onWindowStart:min(eventTon+25,N));
                            [~, onLr] = app.findOneChangepoint(onWindowData); % Likelihood for each point in window.
                            [~, windowTon] = max(onLr); % Snap to most likely changepoint within window.
                            Ton = onWindowTime(windowTon); % Time of new first changepoint.
                            offWindowStart = max(eventToff-25,1); % Start of window for Toff, over which likelihood will be calculated.
                            offWindowData = eventData(offWindowStart:min(eventToff+25,N));
                            offWindowTime = eventTime(offWindowStart:min(eventToff+25,N));
                            [~, offLr] = app.findOneChangepoint(offWindowData); % Likelihood for each point in window.
                            [~, windowToff] = max(offLr); % Snap to most likely changepoint within window.
                            Toff = offWindowTime(windowToff); % Time of new second changepoint.
                        end
                        
                        [~, LrTon] = min(abs(X-Ton)); % First changepoint as index in X, Y, Z (closest point, really, given lower resolution of likelihood).
                        [~, LrToff] = min(abs(X-Toff)); % Second changepoint as index in X, Y, Z (closest point).
                        [~, eventTon] = min(abs(eventTime-Ton)); % First changepoint as index in eventTime.
                        [~, eventToff] = min(abs(eventTime-Toff)); % Second changepoint as index in eventTime.
                        
                        % Calculate the means and standard deviations of the event window (1)
                        % before the first changepoint, (2) between the two changepoints, and
                        % (3) after the second changepoint.
                        y = cumsum(eventData);
                        y2 = cumsum(eventData.^2);
                        x2 = (y(eventToff)-y(eventTon))/(eventToff-eventTon);
                        x13 = (y(N)-y(eventToff)+y(eventTon))/(N-eventToff+eventTon);
                        s213 = ((y2(N)-y2(eventToff)+y2(eventTon)) - ((y(N)-y(eventToff)+y(eventTon)).^2./(N-eventToff+eventTon)))./(N-eventToff+eventTon);
                        s213(s213 < 0) = 0;
                        s22 = ((y2(eventToff)-y2(eventTon))-((y(eventToff)-y(eventTon))^2/(eventToff-eventTon)))/(eventToff-eventTon);
                        s22(s22 < 0) = 0;
                        
                        % Show red lines on the likelihood axes to mark new changepoint
                        % locations.
                        app.changeData(app.analyze.individual.XonPlot, X, Y(LrTon)*ones(M,1), Z(LrTon,:));
                        app.changeData(app.analyze.individual.XoffPlot, X, Y(LrToff)*ones(M,1), Z(LrToff,:));
                        app.changeData(app.analyze.individual.YonPlot, X(LrTon)*ones(M,1), Y, Z(LrTon,:));
                        app.changeData(app.analyze.individual.YoffPlot, X(LrToff)*ones(M,1), Y, Z(LrToff,:));
                        
                        % Set the figure button down callback to pickChngptClick(), and pass store.
                        app.misc.fig.WindowButtonDownFcn = {@app.pickChngptClick store backup};
                        
                        % Show red lines on the individual event axes to mark new changepoint
                        % locations.
                        limy = ylim(app.analyze.individual.eventAxes);
                        app.changeData(app.analyze.individual.TonPlot, [Ton Ton], limy);
                        app.changeData(app.analyze.individual.ToffPlot, [Toff Toff], limy);
                        uistack(app.analyze.individual.TonPlot, 'top');
                        uistack(app.analyze.individual.ToffPlot, 'top');
                        app.changeData(app.analyze.individual.mean1Plot, eventTime(1:eventTon), x13*ones(eventTon,1));
                        app.changeData(app.analyze.individual.upLim1Plot, eventTime(1:eventTon), (x13+sqrt(s213))*ones(eventTon,1));
                        app.changeData(app.analyze.individual.lowLim1Plot, eventTime(1:eventTon), (x13-sqrt(s213))*ones(eventTon,1));
                        app.changeData(app.analyze.individual.mean2Plot, eventTime(eventTon:eventToff), x2*ones(eventToff-eventTon+1,1));
                        app.changeData(app.analyze.individual.upLim2Plot, eventTime(eventTon:eventToff), (x2+sqrt(s22))*ones(eventToff-eventTon+1,1));
                        app.changeData(app.analyze.individual.lowLim2Plot, eventTime(eventTon:eventToff), (x2-sqrt(s22))*ones(eventToff-eventTon+1,1));
                        app.changeData(app.analyze.individual.mean3Plot, eventTime(eventToff:end), x13*ones(N-eventToff+1,1));
                        app.changeData(app.analyze.individual.upLim3Plot, eventTime(eventToff:end), (x13+sqrt(s213))*ones(N-eventToff+1,1));
                        app.changeData(app.analyze.individual.lowLim3Plot, eventTime(eventToff:end), (x13-sqrt(s213))*ones(N-eventToff+1,1));
                    end
                elseif in(2)
                    % If mouse is in the individual event axes:
                    
                    % Get current event's window and likelihood.
                    limy = ylim(app.analyze.individual.eventAxes);
                    
                    % Clicks come in pairs, with each click defining one changepoint.
                    if app.analyze.individual.first
                        % If the next click would mark the first click in a pair:
                        
                        % Assume the click marks the event start changepoint.
                        Ton = p1(2,1); % Time of user-defined first changepoint.
                        if app.analyze.individual.snap.Value
                            % If user wants to use snap to changepoint feature:
                            [~, eventTon] = min(abs(eventTime-Ton)); % User-defined first changepoint as index in eventTime.
                            windowStart = max(eventTon-25,1); % Start of window, over which likelihood will be calculated.
                            windowData = eventData(windowStart:min(eventTon+25,N));
                            windowTime = eventTime(windowStart:min(eventTon+25,N));
                            [~, Lr] = app.findOneChangepoint(windowData); % Likelihood for each point in window.
                            [~, windowTon] = max(Lr); % Snap to most likely changepoint within window.
                            Ton = windowTime(windowTon); % Time of new first changepoint.
                        end
                        [~, LrTon] = min(abs(X-Ton)); % First changepoint as index in X, Y, Z (closest point, really, given lower resolution of likelihood).
                        app.analyze.individual.T1 = Ton;
                        
                        % Show red lines on the individual event and likelihood axes at this
                        % new changepoint.
                        app.changeData(app.analyze.individual.TonPlot, [Ton Ton], limy);
                        app.changeData(app.analyze.individual.XonPlot, X, Y(LrTon)*ones(M,1), Z(LrTon,:));
                        app.changeData(app.analyze.individual.YonPlot, X(LrTon)*ones(M,1), Y, Z(LrTon,:));
                        
                        % Set the figure button down callback to pickChngptClick(), and pass store.
                        app.misc.fig.WindowButtonDownFcn = {@app.pickChngptClick store backup};
                        
                        % Hide all other red lines.
                        app.analyze.individual.XoffPlot.Visible = 'off';
                        app.analyze.individual.YoffPlot.Visible = 'off';
                        app.analyze.individual.ToffPlot.Visible = 'off';
                        app.analyze.individual.mean1Plot.Visible = 'off';
                        app.analyze.individual.upLim1Plot.Visible = 'off';
                        app.analyze.individual.lowLim1Plot.Visible = 'off';
                        app.analyze.individual.mean2Plot.Visible = 'off';
                        app.analyze.individual.upLim2Plot.Visible = 'off';
                        app.analyze.individual.lowLim2Plot.Visible = 'off';
                        app.analyze.individual.mean3Plot.Visible = 'off';
                        app.analyze.individual.upLim3Plot.Visible = 'off';
                        app.analyze.individual.lowLim3Plot.Visible = 'off';
                    else
                        % If the next click would mark the second click in a pair:
                        
                        % First assume the click marks the event end changepoint.
                        Toff = p1(2,1);
                        if app.analyze.individual.snap.Value
                            % If user wants to use snap to changepoint feature:
                            [~, eventToff] = min(abs(eventTime-Toff)); % User-defined second changepoint as index in eventTime.
                            windowStart = max(eventToff-25,1); % Start of window, over which likelihood will be calculated.
                            windowData = eventData(windowStart:min(eventToff+25,N));
                            windowTime = eventTime(windowStart:min(eventToff+25,N));
                            [~, Lr] = app.findOneChangepoint(windowData); % Likelihood for each point in window.
                            [~, windowToff] = max(Lr); % Snap to most likely changepoint within window.
                            Toff = windowTime(windowToff); % Time of new second changepoint.
                        end
                        
                        % Then swap the two changepoints if needed, so that the start
                        % changepoint is less than the end changepoint.
                        Ton = min([app.analyze.individual.T1 Toff]); % Time of first changepoint.
                        Toff = max([app.analyze.individual.T1 Toff]); % Time of second changepoint.
                        [~, LrTon] = min(abs(X-Ton)); % First changepoint as index in X, Y, Z (closest point, really, given lower resolution of likelihood).
                        [~, LrToff] = min(abs(X-Toff)); % Second changepoint as index in X, Y, Z (closest point).
                        [~, eventTon] = min(abs(eventTime-Ton)); % First changepoint as index in eventTime.
                        [~, eventToff] = min(abs(eventTime-Toff)); % Second changepoint as index in eventTime.
                        
                        % Calculate the means and standard deviations of the event window (1)
                        % before the first changepoint, (2) between the two changepoints, and
                        % (3) after the second changepoint.
                        y = cumsum(eventData);
                        y2 = cumsum(eventData.^2);
                        x2 = (y(eventToff)-y(eventTon))/(eventToff-eventTon);
                        x13 = (y(N)-y(eventToff)+y(eventTon))/(N-eventToff+eventTon);
                        s213 = ((y2(N)-y2(eventToff)+y2(eventTon)) - ((y(N)-y(eventToff)+y(eventTon)).^2./(N-eventToff+eventTon)))./(N-eventToff+eventTon);
                        s213(s213 < 0) = 0;
                        s22 = ((y2(eventToff)-y2(eventTon))-((y(eventToff)-y(eventTon))^2/(eventToff-eventTon)))/(eventToff-eventTon);
                        s22(s22 < 0) = 0;
                        
                        % Show red lines on the likelihood axes to mark new changepoint
                        % locations.
                        app.changeData(app.analyze.individual.XonPlot, X, Y(LrTon)*ones(M,1), Z(LrTon,:));
                        app.changeData(app.analyze.individual.XoffPlot, X, Y(LrToff)*ones(M,1), Z(LrToff,:));
                        app.changeData(app.analyze.individual.YonPlot, X(LrTon)*ones(M,1), Y, Z(LrTon,:));
                        app.changeData(app.analyze.individual.YoffPlot, X(LrToff)*ones(M,1), Y, Z(LrToff,:));
                        
                        % Set the figure button down callback to pickChngptClick(), and pass store.
                        app.misc.fig.WindowButtonDownFcn = {@app.pickChngptClick store backup};

                        % Show red lines on the individual event axes to mark new changepoint
                        % locations.
                        app.changeData(app.analyze.individual.TonPlot, [Ton Ton], limy);
                        app.changeData(app.analyze.individual.ToffPlot, [Toff Toff], limy);
                        app.changeData(app.analyze.individual.mean1Plot, eventTime(1:eventTon), x13*ones(eventTon,1));
                        app.changeData(app.analyze.individual.upLim1Plot, eventTime(1:eventTon), (x13+sqrt(s213))*ones(eventTon,1));
                        app.changeData(app.analyze.individual.lowLim1Plot, eventTime(1:eventTon), (x13-sqrt(s213))*ones(eventTon,1));
                        app.changeData(app.analyze.individual.mean2Plot, eventTime(eventTon:eventToff), x2*ones(eventToff-eventTon+1,1));
                        app.changeData(app.analyze.individual.upLim2Plot, eventTime(eventTon:eventToff), (x2+sqrt(s22))*ones(eventToff-eventTon+1,1));
                        app.changeData(app.analyze.individual.lowLim2Plot, eventTime(eventTon:eventToff), (x2-sqrt(s22))*ones(eventToff-eventTon+1,1));
                        app.changeData(app.analyze.individual.mean3Plot, eventTime(eventToff:end), x13*ones(N-eventToff+1,1));
                        app.changeData(app.analyze.individual.upLim3Plot, eventTime(eventToff:end), (x13+sqrt(s213))*ones(N-eventToff+1,1));
                        app.changeData(app.analyze.individual.lowLim3Plot, eventTime(eventToff:end), (x13-sqrt(s213))*ones(N-eventToff+1,1));
                    end
                else
                    % Disable the figure button down callback.
                    app.misc.fig.WindowButtonDownFcn = [];

                    % If mouse is in neither axes, hide all red lines.
                    app.analyze.individual.XonPlot.Visible = 'off';
                    app.analyze.individual.XoffPlot.Visible = 'off';
                    app.analyze.individual.YonPlot.Visible = 'off';
                    app.analyze.individual.YoffPlot.Visible = 'off';
                    app.analyze.individual.TonPlot.Visible = 'off';
                    app.analyze.individual.ToffPlot.Visible = 'off';
                    app.analyze.individual.mean1Plot.Visible = 'off';
                    app.analyze.individual.upLim1Plot.Visible = 'off';
                    app.analyze.individual.lowLim1Plot.Visible = 'off';
                    app.analyze.individual.mean2Plot.Visible = 'off';
                    app.analyze.individual.upLim2Plot.Visible = 'off';
                    app.analyze.individual.lowLim2Plot.Visible = 'off';
                    app.analyze.individual.mean3Plot.Visible = 'off';
                    app.analyze.individual.upLim3Plot.Visible = 'off';
                    app.analyze.individual.lowLim3Plot.Visible = 'off';
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                if isequal(app.misc.fig.WindowButtonMotionFcn, thisMotionCallback)
                    app.misc.fig.WindowButtonMotionFcn = [];
                    app.unexpectedError(ME);
                    app.restoreBackupTab1(backup);
                end
            end
        end
        
        function pickChngptClick(app, ~, ~, store, backup)
            % Button down callback for the figure when user is setting manual changepoints
            % by clicking the likelihood plot or individual event plot.
            
            try
                if ~isequal(app.misc.tab1, app.misc.tabgroup.SelectedTab)
                    return
                end

                in = app.inAxes([app.analyze.individual.likelihoodAxes...
                    app.analyze.individual.eventAxes]); % Whether in either axes.
                if in(1) || (in(2) && ~app.analyze.individual.first)
                    % If in likelihood axes OR in individual event axes and on the second
                    % click:

                    app.misc.fig.WindowButtonMotionFcn = app.misc.storeMotionCallback; % Reset the figure motion callback.
                    app.misc.fig.WindowButtonDownFcn = app.misc.storeClickCallback; % Reset the figure button down callback.
                    app.misc.fig.WindowKeyPressFcn = @app.nullKeyPressCallback;

                    % Update plots on likelihood axes.
                    app.changeData(app.analyze.individual.defXonPlot, app.analyze.individual.XonPlot);
                    app.changeData(app.analyze.individual.defXoffPlot, app.analyze.individual.XoffPlot);
                    app.changeData(app.analyze.individual.defYonPlot, app.analyze.individual.YonPlot);
                    app.changeData(app.analyze.individual.defYoffPlot, app.analyze.individual.YoffPlot); 

                    % Update plots on individual event axes.
                    app.changeData(app.analyze.individual.defTonPlot, app.analyze.individual.TonPlot);
                    app.changeData(app.analyze.individual.defToffPlot, app.analyze.individual.ToffPlot);
                    app.changeData(app.analyze.individual.defMean1Plot, app.analyze.individual.mean1Plot);
                    app.changeData(app.analyze.individual.defUpLim1Plot, app.analyze.individual.upLim1Plot);
                    app.changeData(app.analyze.individual.defLowLim1Plot, app.analyze.individual.lowLim1Plot);
                    app.changeData(app.analyze.individual.defMean2Plot, app.analyze.individual.mean2Plot);
                    app.changeData(app.analyze.individual.defUpLim2Plot, app.analyze.individual.upLim2Plot);
                    app.changeData(app.analyze.individual.defLowLim2Plot, app.analyze.individual.lowLim2Plot);
                    app.changeData(app.analyze.individual.defMean3Plot, app.analyze.individual.mean3Plot);
                    app.changeData(app.analyze.individual.defUpLim3Plot, app.analyze.individual.upLim3Plot);
                    app.changeData(app.analyze.individual.defLowLim3Plot, app.analyze.individual.lowLim3Plot);
                    
                    % Update allEvents to reflect new changepoints.
                    j = app.getCurrentEventAsIndex;
                    fileStart = app.analyze.individual.eventAxes.XLim(1);
                    fileStop = app.analyze.individual.eventAxes.XLim(2);
                    [~, app.data.windowStarts(j)] = min(abs(app.data.time-fileStart));
                    [~, app.data.windowStops(j)] = min(abs(app.data.time-fileStop));
                    app.resetAxesDefaultLimits(app.analyze.individual.eventAxes);
                    app.allEvents(j,:) =...
                        [find(app.data.time==app.analyze.individual.XonPlot.YData(1))...
                        find(app.data.time==app.analyze.individual.XoffPlot.YData(1))];
                    if ~isempty(app.misc.listenerError)
                        rethrow(app.misc.listenerError)
                    end
                    
                    app.analyze.individual.likelihoodAxes.XColor = 'k';
                    app.analyze.individual.likelihoodAxes.YColor = 'k';
                    app.analyze.individual.likelihoodAxes.ZColor = 'k';
                    app.analyze.individual.eventAxes.XColor = 'k';
                    app.analyze.individual.eventAxes.YColor = 'k';
                    app.analyze.individual.eventAxes.ZColor = 'k';
                    
                    app.restore(store);
                elseif in(2)
                    % If in individual event axes and on the first click, indicate that user
                    % is now on the second click.
                    app.analyze.individual.first = false;
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowButtonDownFcn = [];
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function pickChngptEscape(app, ~, event, store, backup)
            % When pick changepoint mode is active, this is called whenever user presses
            % any key. If user presses the escape key, program will exit pick changepoint
            % mode and return to state prior to entering the mode.
            
            try
                if strcmp(event.Key, 'escape')
                    app.misc.fig.WindowButtonMotionFcn = app.misc.storeMotionCallback; % Reset the figure motion callback.
                    app.misc.fig.WindowButtonDownFcn = app.misc.storeClickCallback; % Reset the figure button down callback.
                    app.misc.fig.WindowKeyPressFcn = @app.nullKeyPressCallback; % Reset the figure key press callback.

                    % Leave black lines where they are, and simply hide red lines.
                    app.analyze.individual.XonPlot.Visible = 'off';
                    app.analyze.individual.XoffPlot.Visible = 'off';
                    app.analyze.individual.YonPlot.Visible = 'off';
                    app.analyze.individual.YoffPlot.Visible = 'off';
                    app.analyze.individual.TonPlot.Visible = 'off';
                    app.analyze.individual.ToffPlot.Visible = 'off';
                    app.analyze.individual.mean1Plot.Visible = 'off';
                    app.analyze.individual.upLim1Plot.Visible = 'off';
                    app.analyze.individual.lowLim1Plot.Visible = 'off';
                    app.analyze.individual.mean2Plot.Visible = 'off';
                    app.analyze.individual.upLim2Plot.Visible = 'off';
                    app.analyze.individual.lowLim2Plot.Visible = 'off';
                    app.analyze.individual.mean3Plot.Visible = 'off';
                    app.analyze.individual.upLim3Plot.Visible = 'off';
                    app.analyze.individual.lowLim3Plot.Visible = 'off';

                    % Reset axes colors.
                    app.analyze.individual.likelihoodAxes.XColor = 'k';
                    app.analyze.individual.likelihoodAxes.YColor = 'k';
                    app.analyze.individual.likelihoodAxes.ZColor = 'k';
                    app.analyze.individual.eventAxes.XColor = 'k';
                    app.analyze.individual.eventAxes.YColor = 'k';
                    app.analyze.individual.eventAxes.ZColor = 'k';
                    
                    app.restore(store);
                end
            catch ME
                if ~isvalid(app)
                    return
                end
                app.misc.fig.WindowKeyPressFcn = [];
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function currentEventCallback(app, src, ~)
            % Called when user changes the value of currentEventEdit. Makes sure input is
            % valid.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get input.
                d = str2double(src.String);

                % Get current value.
                prev = int2str(app.analyze.individual.currentEvent);

                if ~isfinite(d) || ~isreal(d)
                    % If input is infinite or NaN, reset edit text string to current value and
                    % return.
                    uiwait(errordlg('The current event must be an integer.', 'Error Dialog', 'modal'))
                    src.String = prev;
                    app.restore(store);
                    return
                end

                if mod(d,1) ~= 0
                    % If input is a decimal, round down to nearest integer.
                    warndlg('The current event must be an integer. Rounding down.')
                    d = floor(d);
                end

                switch app.showRemoved
                    case 0
                        lim = sum(app.activeEvents);
                    case 1
                        lim = length(app.activeEvents);
                end
                if d < 1 || d > lim
                    % If input is not within range, reset edit text string to current value
                    % and return.
                    uiwait(errordlg(sprintf('The current event should be within %d and %d.', 1, lim), 'Error Dialog', 'modal'))
                    src.String = prev;
                    app.restore(store);
                    return
                end

                % Input is an integer within range here.
                if app.changeCurrentEvent(d)
                    % Change current event to new value, and if new value does not match
                    % previous value, update plots.
                    app.plotIndividualEvent;
                    app.updateEventVisNColor;
                end
                
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.String = int2str(app.analyze.individual.currentEvent);
            end
        end
        
        function j = getCurrentEventAsIndex(app)
            % Returns the location of app.analyze.individual.currentEvent in
            % allEvents/activeEvents
            %
            % Suppose activeEvents = [1 1 1 0 1 0 1] and
            % app.analyze.individual.currentEvent = 4.
            %
            % Event numbering depends on showRemoved:
            %
            %     if showRemoved:  1 2 3 4 5 6 7
            %       activeEvents: [1 1 1 0 1 0 1]
            % if not showRemoved:  1 2 3   4   5
            %
            % If showRemoved is true, the current event is the fourth event, but if
            % showRemoved is false, the current event is the fifth event.
            
            i = app.analyze.individual.currentEvent;
            switch app.showRemoved
                case 0
                    % If showRemoved is false, the index of the current event might not be
                    % equal to currentEvent.
                    j = find(app.activeEvents, i);
                    if ~isempty(j)
                       j = j(end);
                    end
                case 1
                    % If showRemoved is true, nothing changes.
                    j = i;
            end
        end
        
        function removeEventCallback(app, src, ~)
            % Called when user clicks individual.removeBtn. Updates activeEvents
            % accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                i = app.getCurrentEventAsIndex;
                v = src.Value;
                app.analyze.individual.removedByUser(i) = v;
                app.analyze.individual.includedByUser(i) = ~v;
                app.analyze.individual.storeCurrentEvent = app.getCurrentEventAsIndex;
                app.activeEvents(i) = ~v;
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.Value = ~app.activeEvents(i);
            end
        end
        
        function nextEventCallback(app, ~, ~)
            % Called when user clicks nextEventBtn. Increases the current event by 1 and
            % updates plots, if needed.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                if app.changeCurrentEvent(app.analyze.individual.currentEvent+1)
                    % Change current event to new value, and if new value does not match
                    % previous value, update plots.
                    app.plotIndividualEvent;
                    app.updateEventVisNColor;
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function prevEventCallback(app, ~, ~)
            % Called when user clicks nextEventBtn. Decreases the current event by 1 and
            % updates plots, if needed.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                if app.changeCurrentEvent(app.analyze.individual.currentEvent-1)
                    % Change current event to new value, and if new value does not match
                    % previous value, update plots.
                    app.plotIndividualEvent;
                    app.updateEventVisNColor;
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function showRemovedCallback(app, src, ~)
            % Called when user clicks showRemovedBtn. Updates showRemoved accordingly.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;
                app.analyze.individual.storeCurrentEvent = app.getCurrentEventAsIndex;
                app.showRemoved = src.Value;
                if ~isempty(app.misc.listenerError)
                    rethrow(app.misc.listenerError)
                end
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
                src.Value = app.showRemoved;
            end
        end
        
        function needUpdate = changeCurrentEvent(app, target)
            % Sets app.analyze.individual.currentEvent to target and returns whether plots
            % need to be updated.
            
            from = app.getCurrentEventAsIndex; % Get current event before app.analyze.individual.currentEvent changes.
            
            app.analyze.individual.currentEvent = target; % Update currentEvent.
            app.analyze.individual.currentEventEdit.String = int2str(target); % Update the edit text.
            
            app.updateCheckboxAndLabel;
            
            % Determine whether prevEventBtn and nextEventBtn should be enabled.
            switch app.showRemoved
                case 0
                    lim = sum(app.activeEvents);
                case 1
                    lim = length(app.activeEvents);
            end
            app.addToNew(target > 1, app.analyze.individual.prevEventBtn);
            app.addToNew(target < lim, app.analyze.individual.nextEventBtn);
            
            % Set needUpdate to reflect whether the new current event is different from
            % the previous current event and thus whether plots need updating.
            needUpdate = app.getCurrentEventAsIndex ~= from;
            % This value is only useful if activeEvents and showRemoved have not changed
            % since the last time currentEvent was set. postSetCallback() handles the
            % cases where activeEvents or showRemoved do change.
        end
        
        function updateCheckboxAndLabel(app)
            % Updates the remove event checkbox and label indicating the reason for
            % removal to reflect activeEvents, currentEvent, and removedBy arrays.
            
            ce = app.getCurrentEventAsIndex;
            app.analyze.individual.removeBtn.Value = ~app.activeEvents(ce); % Update the checkbox.
            
            % Update the label.
            s = '';
            if app.analyze.individual.removedByDur(ce), s = [s 'dur + ']; end
            if app.analyze.individual.removedBySep(ce), s = [s 'sep + ']; end
            if app.analyze.individual.removedByUser(ce), s = [s 'user']; end
            s = strip(strip(strip(s), '+'));
            if ~app.activeEvents(ce)
                app.analyze.individual.label.String = ['Removed: ' s];
                app.analyze.individual.label.ForegroundColor = app.misc.colors.removeLabel;
            elseif app.analyze.individual.includedByUser(ce) && (app.analyze.individual.removedByDur(ce)...
                    || app.analyze.individual.removedBySep(ce))
                app.analyze.individual.label.String = 'Included by user';
                app.analyze.individual.label.ForegroundColor = app.misc.colors.includeLabel;
            else
                app.analyze.individual.label.String = '';
            end
            
            drawnow;
        end
        
        function createExcelCallback(app, ~, ~)
            % Called when user clicks createBtn. Allows user to create a blank excel
            % document and add it to the excel listbox.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Let user specify name for new file.
                [name, path] = uiputfile('sample.xlsx');
                if name == 0
                    app.restore(store);
                    return
                end
                [~, name] = fileparts([path name]);
                name = [name '.xlsx'];
                app.updateCurrentTask(sprintf('Creating %s.', fullfile(path, name)));

                % Add new file's name and path to the list of excel files.
                if app.excel.files(1,1) == ""
                    app.excel.files(1,[1 2]) = [string(path) string(name)];
                else
                    app.excel.files(end+1,[1 2]) = [string(path) string(name)];
                end
                app.excel.listbox.String = app.excel.files(:,2); % Update the listbox.

                % File is not actually created up to this point. This line creates the file.
                writetable(table(""), fullfile(path, name), 'WriteVariableNames', 0);

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function findExcelCallback(app, ~, ~)
            % Called when user clicks findBtn. Allows user to find an existing excel
            % document and add it to the excel listbox.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Let user pick .xlsx file.
                [name, path] = uigetfile('.xlsx');
                if name == 0
                    app.restore(store);
                    return
                end
                [~, ~, ext] = fileparts([path name]);
                if ~strcmp(ext, '.xlsx')
                    uiwait(errordlg('You must choose a .xlsx file.', 'Error Dialog', 'modal'))
                    app.restore(store);
                    return
                end
                app.updateCurrentTask(sprintf('Adding %s to the list.', fullfile(path, name)));

                % Add this file's name and path to the list of excel files.
                if app.excel.files(1,1) == ""
                    app.excel.files(1,[1 2]) = [string(path) string(name)];
                else
                    if any(ismember(app.excel.files, [string(path) string(name)], 'rows'))
                        % If chosen file is already in list, exit.
                        app.restore(store);
                        return
                    end
                    app.excel.files(end+1,[1 2]) = [string(path) string(name)];
                end
                app.excel.listbox.String = app.excel.files(:,2); % Update the listbox.

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function addToExcelCallback(app, ~, ~)
            % Called when user clicks addBtn. Writes pertinent information about current
            % .txt file and results from analysis to all selected excel files.
            
            backup = app.storeBackupTab1;
            try
                store = app.disable;
                drawnow;

                % Get selected excel files.
                v = app.excel.listbox.Value; % Indices of selected excel files.
                if isempty(v)
                    % If no excel files are selected, exit.
                    app.restore(store);
                    return
                end
                app.updateCurrentTask('Preparing to save data to selected .xlsx files.');

                paths = app.excel.files(v,1); % Paths of selected excel files
                files = app.excel.files(v,2); % Names of selected excel files.

                % Determine data to write to file.
                indices = sort([find(diff([0; ~app.controls.deselect.indices])>0);...
                    find(diff([~app.controls.deselect.indices; 0])<0)]); % Start and end indices of each deselect patch.
                if isequal(app.misc.active.allEvents.whichBeads, app.controls.beadABtn)
                    % Just bead A.
                    beads = 'only A';
                elseif isequal(app.misc.active.allEvents.whichBeads, app.controls.beadBBtn)
                    % Just bead B.
                    beads = 'only B';
                else
                    % Both.
                    beads = 'both';
                end
                sep = app.getActiveSep; % Event separations.
                toFile = [app.getActiveEvents...        % Indices before event start and before event end.
                    app.getActiveDur...                 % Event durations in points.
                    sep(1:end-1)...                     % Preceding separations in points.
                    sep(2:end)...                       % Succeeding separations in points.
                    app.getActiveEvents/app.data.fs...  % Times of event start and end times in seconds.
                    app.getActiveDur/app.data.fs...     % Event durations in seconds.
                    sep(1:end-1)/app.data.fs...         % Preceding separations in seconds.
                    sep(2:end)/app.data.fs...           % Succeeding separations in seconds.
                    app.prepExcel];                     % Output from prepExcel().
                colNames = {'index before start (# pts)'...
                    'index before end (# pts)'...
                    'duration (# pts)'...
                    'before (# pts)'...
                    'after (# pts)'...
                    'time before start (s)'...
                    'time before end (s)'...
                    'duration (s)'...
                    'before (s)'...
                    'after (s)'...
                    'A pos before attach (nm)'...
                    'A pos after attach (nm)'...
                    'A pos before detach (nm)'...
                    'A pos after detach (nm)'...
                    'A step 1 (nm)'...
                    'A step 2 (nm)'...
                    'A total step (nm)'...
                    'B pos before attach (nm)'...
                    'B pos after attach (nm)'...
                    'B pos before detach (nm)'...
                    'B pos after detach (nm)'...
                    'B step 1 (nm)'...
                    'B step 2 (nm)'...
                    'B total step (nm)'...
                    'A avg force before (pN)'...
                    'A avg force during (pN)'...
                    'A force before detach (pN)'...
                    'A force after detach (pN)'...
                    'A avg force after (pN)'...
                    'B avg force before (pN)'...
                    'B avg force during (pN)'...
                    'B force before detach (pN)'...
                    'B force after detach (pN)'...
                    'B avg force after (pN)'};

                notFound = false(size(paths, 1),1);
                for i = 1:size(paths, 1)
                    % For each excel file...
                    fullpath = [char(paths(i)) char(files(i))];
                    if ~isfile(fullpath)
                        notFound(i) = true;
                        % ...if file not found, skip it.
                        continue;
                    end
                    app.updateCurrentTask(sprintf('Saving data to %s.', fullpath));
                    s = app.getNextSheet(fullpath); % ...determine which sheet to write to...
                    fileID = fopen(fullpath, 'r+'); % ...try to open the file for writing...
                    if fileID ~= -1
                        % ...if able to open the file for writing...
                        fclose(fileID); % ...close it, because writetable opens it independently...

                        % ...write data to the file, with sheet determined by s, column
                        % determined by getExcelCol(), and row determined by values stored in
                        % app.excel.cell...
                        warning('off', 'MATLAB:xlswrite:AddSheet')
                        writetable(table({strcat(app.load.path, '/', app.load.name)}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.path], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'sampling frequency (Hz)', app.data.fs}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.fs], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'number of events', size(toFile, 1)}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.N], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'flipped', app.flipped}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.flipped], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'ignore data between these indices (start 1, stop 1, start 2, stop 2, etc)', indices}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.deselectIndices], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'minimum duration', app.minDur}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.minDur], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'minimum separation', app.minSep}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.minSep], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'moving average filter window', app.controls.windows.covwindow}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.covwindow], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'2nd order Savitzky-Golay filter window', app.controls.windows.covsmooth}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.covsmooth], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'auto peak 1', app.autop1}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.autopeak1], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'auto peak 2', app.autop2}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.autopeak2], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'auto min', app.automin}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.automin], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'manual peak 1', app.manp1}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.manpeak1], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'manual peak 2', app.manp2}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.manpeak2], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'manual min', app.manmin}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.manmin], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'use minimum', app.usemin}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.usemin], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'use manual', app.useman}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.useman], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'use change point', logical(app.misc.active.allEvents.useChangepoint)}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.useChngpt], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table({'analyzed bead(s)', beads}), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.whichBeads], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table(colNames), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.colNames], 'Sheet', s, 'WriteVariableNames', 0);
                        writetable(table(toFile), fullpath, 'Range', [app.getExcelCol(1) app.excel.cell.data], 'Sheet', s, 'WriteVariableNames', 0);
                        warning('on', 'MATLAB:xlswrite:AddSheet')
                    else
                        % ...if unable to open the file for writing, skip it.
                        notFound(i) = true;
                    end
                end
                app.updateCurrentTask('Cleaning up.');
                app.excel.listbox.Value = [];
                if any(notFound)
                    warndlg(['The following files were not found or could not be opened:'...
                        sprintf('\n    %s', files{notFound})...
                        sprintf('\nThey have been removed from the list.')])
                    app.excel.files(notFound,:) = [];
                end
                if size(app.excel.files, 1) == 0
                    app.excel.files = strings;
                    app.excel.listbox.String = strings;
                else
                    app.excel.listbox.String = app.excel.files(:,2); % Update the listbox.
                end
                app.restore(store);
            catch ME
                try
                    fclose(fileID);
                catch
                end
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab1(backup);
            end
        end
        
        function s = getNextSheet(app, fullpath)
            % Determines the next blank sheet in the Excel workbook specified by fullpath.
            
            nSheets = app.numSheets(fullpath);
            for i = 1:nSheets
                % For each sheet...
                T = readtable(fullpath, 'Sheet', i, 'Range', 'A1:A1', 'ReadVariableNames', 0);
                if ~iscell(T{:,:}) && isnan(T{:,:})
                    % ...if cell A1 is empty, this sheet should be written to next. Exit.
                    s = i;
                    return
                end
            end
            % If no existing sheets had an empty A1 cell, write to a new sheet.
            s = nSheets+1;
        end
        
        function toExcel = prepExcel(app)
            % Returns data to write to excel files. For each event, returns:
            %  Position of bead A before and after binding and unbinding.
            %  Position of bead B before and after binding and unbinding.
            %  Bead A step 1, step 2, and total step.
            %  Bead B step 1, step 2, and total step.
            %  Force on bead A just before, during, and just after event.
            %  Force on bead B just before, during, and just after event.
            
            [posA, stepA] = app.getPosNStepSizes('A');
            [posB, stepB] = app.getPosNStepSizes('B');
            forceA = app.getForce('A');
            forceB = app.getForce('B');
            toExcel = [posA...
                stepA...
                posB...
                stepB...
                forceA...
                forceB];
        end
    end
    
    % Private tab 2 methods
    methods (Access = private)
        
        function backup = storeBackupTab2(app)
            % Copies by value - not by reference - all tab 2 properties which could
            % potentially change.
            
            backup.misc.fig.Position = app.misc.fig.Position;
            
            backup.load2.inputBtn = storeButton(app.load2.inputBtn);
            
            backup.load2.table.Units = app.load2.table.Units;
            backup.load2.table.Position = app.load2.table.Position;
            backup.load2.tableData = app.load2.tableData;
            backup.load2.table.Visible = app.load2.table.Visible;
            backup.load2.table.Enable = app.load2.table.Enable;
            
            backup.data2.fs = app.data2.fs;
            backup.data2.events = app.data2.events;
            backup.data2.beadA = app.data2.beadA;
            backup.data2.beadB = app.data2.beadB;
            backup.load2.beadABtn = storeButton(app.load2.beadABtn);
            backup.load2.beadBBtn = storeButton(app.load2.beadBBtn);
            backup.load2.bothBeadsBtn = storeButton(app.load2.bothBeadsBtn);
            backup.load2.whichBeadBtnGrp.SelectedObject = app.load2.whichBeadBtnGrp.SelectedObject;
            backup.load2.whichBeadBtnGrp.Visible = app.load2.whichBeadBtnGrp.Visible;
            backup.load2.skipBtn = storeButton(app.load2.skipBtn);
            
            backup.analyze2.misc.miscAxes = storeAxes(app.analyze2.misc.miscAxes);
            backup.analyze2.misc.miscAxes.YScale = app.analyze2.misc.miscAxes.YScale;
            backup.analyze2.misc.miscAxes.YLabel.String = app.analyze2.misc.miscAxes.YLabel.String;
            backup.analyze2.misc.miscAxes.XLabel.String = app.analyze2.misc.miscAxes.XLabel.String;
            backup.analyze2.misc.miscPlot = storePlot(app.analyze2.misc.miscPlot);
            backup.analyze2.misc.miscPlot.Marker = app.analyze2.misc.miscPlot.Marker;
            backup.analyze2.misc.miscPlot.LineStyle = app.analyze2.misc.miscPlot.LineStyle;
            backup.analyze2.misc.miscFitPlot = storePlot(app.analyze2.misc.miscFitPlot);
            backup.analyze2.misc.miscLabel = storeText(app.analyze2.misc.miscLabel);
            backup.analyze2.misc.data = app.analyze2.misc.data;
            backup.analyze2.misc.kDur = app.analyze2.misc.kDur;
            backup.analyze2.misc.meanStep1 = app.analyze2.misc.meanStep1;
            backup.analyze2.misc.stdStep1 = app.analyze2.misc.stdStep1;
            backup.analyze2.misc.meanStep2 = app.analyze2.misc.meanStep2;
            backup.analyze2.misc.stdStep2 = app.analyze2.misc.stdStep2;
            backup.analyze2.misc.meanTotalStep = app.analyze2.misc.meanTotalStep;
            backup.analyze2.misc.stdTotalStep = app.analyze2.misc.stdTotalStep;
            
            backup.data2.durs = app.data2.durs;
            backup.data2.beadAStep1 = app.data2.beadAStep1;
            backup.data2.beadBStep1 = app.data2.beadBStep1;
            backup.data2.beadAStep2 = app.data2.beadAStep2;
            backup.data2.beadBStep2 = app.data2.beadBStep2;
            backup.data2.beadATotalStep = app.data2.beadATotalStep;
            backup.data2.beadBTotalStep = app.data2.beadBTotalStep;
            backup.data2.forceA = app.data2.forceA;
            backup.data2.forceB = app.data2.forceB;
            
            backup.analyze2.ensemble.ensembleAxes = storeAxes(app.analyze2.ensemble.ensembleAxes);
            backup.analyze2.ensemble.ensembleAxes.YLabel.String = app.analyze2.ensemble.ensembleAxes.YLabel.String;
            backup.analyze2.ensemble.globalFit = app.analyze2.ensemble.globalFit;
            backup.analyze2.ensemble.start = app.analyze2.ensemble.start;
            backup.analyze2.ensemble.end = app.analyze2.ensemble.end;
            backup.analyze2.ensemble.startFit = app.analyze2.ensemble.startFit;
            backup.analyze2.ensemble.endFit = app.analyze2.ensemble.endFit;
            backup.analyze2.ensemble.p = app.analyze2.ensemble.p;
            backup.analyze2.ensemble.startPlot = storePlot(app.analyze2.ensemble.startPlot);
            backup.analyze2.ensemble.endPlot = storePlot(app.analyze2.ensemble.endPlot);
            backup.analyze2.ensemble.startFitPlot = storePlot(app.analyze2.ensemble.startFitPlot);
            backup.analyze2.ensemble.endFitPlot = storePlot(app.analyze2.ensemble.endFitPlot);
            backup.analyze2.ensemble.startLabel = storeText(app.analyze2.ensemble.startLabel);
            backup.analyze2.ensemble.endLabel = storeText(app.analyze2.ensemble.endLabel);
            backup.analyze2.ensemble.stepLabel = storeText(app.analyze2.ensemble.stepLabel);
            
            backup.misc.toBeEnabled = app.misc.toBeEnabled;
            backup.misc.disabled = app.misc.disabled;
            backup.misc.skip = app.misc.skip;
            enableable = app.getEnableable;
            for e = 1:numel(enableable)
                backup.enableable(e).Enable = enableable(e).Enable;
            end
            backup.load.currentTask = storeText(app.load.currentTask);
            backup.load2.currentTask = storeText(app.load2.currentTask);
            if ~app.misc.axToolbars
                backup.misc.menus.Enable = app.misc.menus(1).Enable;
            end
            
            function copyB = storeButton(b)
                copyB.Visible = b.Visible;
                copyB.Enable = b.Enable;
            end
            
            function copyT = storeText(t)
                copyT.Visible = t.Visible;
                if isprop(t, 'Enable'), copyT.Enable = t.Enable; end
                copyT.String = t.String;
                copyT.Position = t.Position;
            end
            
            function copyA = storeAxes(a)
                copyA.Visible = a.Visible;
                copyA.XColor = a.XColor;
                copyA.YColor = a.YColor;
                copyA.ZColor = a.ZColor;
                copyA.Position = a.Position;
                copyA.View = a.View;
                copyA.Limits = axis(a);
                if app.misc.axToolbars
                    copyA.Toolbar.Visible = a.Toolbar.Visible;
                end
            end
            
            function copyP = storePlot(p)
                if isempty(p), copyP = []; return, end
                copyP(numel(p)) = struct;
                for i = 1:numel(p)
                    copyP(i).Visible = p(i).Visible;
                    copyP(i).XData = p(i).XData;
                    copyP(i).YData = p(i).YData;
                    if isprop(p(i), 'ZData'), copyP(i).ZData = p(i).ZData; end
                    copyP(i).Parent = p(i).Parent;
                end
            end
        end
        
        function restoreBackupTab2(app, backup, mode)
            % Restores tab 2 properties to the values contained in the struct backup. The
            % subset of properties which are updated depends on mode.
            
            if nargin < 3
                mode = 'all';
            end
            
            generalQ = any(strcmp(mode, {'all'}));
            
            restoreButton(app.load2.inputBtn, backup.load2.inputBtn);
            
            app.load2.table.Units = backup.load2.table.Units;
            app.load2.table.Position = backup.load2.table.Position;
            if strcmp(backup.load2.table.Visible, 'on')
                if ~isempty(backup.load2.tableData)
                    app.load2.tWidths = app.fillTable(app.load2.table, backup.load2.tableData, app.load2.maxHeight);
                else
                    app.load2.table.Data = backup.load2.tableData;
                end
            end
            app.load2.table.Visible = backup.load2.table.Visible;
            app.load2.table.Enable = backup.load2.table.Enable;
              
            app.data2.fs = backup.data2.fs;
            app.data2.events = backup.data2.events;
            app.data2.beadA = backup.data2.beadA;
            app.data2.beadB = backup.data2.beadB;
            restoreButton(app.load2.beadABtn, backup.load2.beadABtn);
            restoreButton(app.load2.beadBBtn, backup.load2.beadBBtn);
            restoreButton(app.load2.bothBeadsBtn, backup.load2.bothBeadsBtn);
            app.load2.whichBeadBtnGrp.SelectedObject = backup.load2.whichBeadBtnGrp.SelectedObject;
            app.load2.whichBeadBtnGrp.Visible = backup.load2.whichBeadBtnGrp.Visible;
            restoreButton(app.load2.skipBtn, backup.load2.skipBtn);
            
            restoreAxes(app.analyze2.misc.miscAxes, backup.analyze2.misc.miscAxes);
            app.analyze2.misc.miscAxes.YScale = backup.analyze2.misc.miscAxes.YScale;
            app.analyze2.misc.miscAxes.YLabel.String = backup.analyze2.misc.miscAxes.YLabel.String;
            app.analyze2.misc.miscAxes.XLabel.String = backup.analyze2.misc.miscAxes.XLabel.String;
            restorePlot(app.analyze2.misc.miscPlot, backup.analyze2.misc.miscPlot);
            app.analyze2.misc.miscPlot.Marker = backup.analyze2.misc.miscPlot.Marker;
            app.analyze2.misc.miscPlot.LineStyle = backup.analyze2.misc.miscPlot.LineStyle;
            restorePlot(app.analyze2.misc.miscFitPlot, backup.analyze2.misc.miscFitPlot);
            restoreText(app.analyze2.misc.miscLabel, backup.analyze2.misc.miscLabel);
            app.analyze2.misc.data = backup.analyze2.misc.data;
            app.analyze2.misc.kDur = backup.analyze2.misc.kDur;
            app.analyze2.misc.meanStep1 = backup.analyze2.misc.meanStep1;
            app.analyze2.misc.stdStep1 = backup.analyze2.misc.stdStep1;
            app.analyze2.misc.meanStep2 = backup.analyze2.misc.meanStep2;
            app.analyze2.misc.stdStep2 = backup.analyze2.misc.stdStep2;
            app.analyze2.misc.meanTotalStep = backup.analyze2.misc.meanTotalStep;
            app.analyze2.misc.stdTotalStep = backup.analyze2.misc.stdTotalStep;
            
            app.data2.durs = backup.data2.durs;
            app.data2.beadAStep1 = backup.data2.beadAStep1;
            app.data2.beadBStep1 = backup.data2.beadBStep1;
            app.data2.beadAStep2 = backup.data2.beadAStep2;
            app.data2.beadBStep2 = backup.data2.beadBStep2;
            app.data2.beadATotalStep = backup.data2.beadATotalStep;
            app.data2.beadBTotalStep = backup.data2.beadBTotalStep;
            app.data2.forceA = backup.data2.forceA;
            app.data2.forceB = backup.data2.forceB;
            
            restoreAxes(app.analyze2.ensemble.ensembleAxes, backup.analyze2.ensemble.ensembleAxes);
            app.analyze2.ensemble.ensembleAxes.YLabel.String = backup.analyze2.ensemble.ensembleAxes.YLabel.String;
            app.analyze2.ensemble.globalFit = backup.analyze2.ensemble.globalFit;
            app.analyze2.ensemble.start = backup.analyze2.ensemble.start;
            app.analyze2.ensemble.end = backup.analyze2.ensemble.end;
            app.analyze2.ensemble.startFit = backup.analyze2.ensemble.startFit;
            app.analyze2.ensemble.endFit = backup.analyze2.ensemble.endFit;
            app.analyze2.ensemble.p = backup.analyze2.ensemble.p;
            restorePlot(app.analyze2.ensemble.startPlot, backup.analyze2.ensemble.startPlot);
            restorePlot(app.analyze2.ensemble.endPlot, backup.analyze2.ensemble.endPlot);
            restorePlot(app.analyze2.ensemble.startFitPlot, backup.analyze2.ensemble.startFitPlot);
            restorePlot(app.analyze2.ensemble.endFitPlot, backup.analyze2.ensemble.endFitPlot);
            restoreText(app.analyze2.ensemble.startLabel, backup.analyze2.ensemble.startLabel);
            restoreText(app.analyze2.ensemble.endLabel, backup.analyze2.ensemble.endLabel);
            restoreText(app.analyze2.ensemble.stepLabel, backup.analyze2.ensemble.stepLabel);
            
            % So that 'Reset/Restore View' buttons work properly, we had to update all
            % plots while XLimMode, YLimMode, etc were 'auto'. After setting the limits,
            % these properties will become 'manual'.
            restoreAxesLimits(app.analyze2.misc.miscAxes, backup.analyze2.misc.miscAxes);
            restoreAxesLimits(app.analyze2.ensemble.ensembleAxes, backup.analyze2.ensemble.ensembleAxes);
            
            if generalQ
                app.misc.fig.Position = backup.misc.fig.Position;
                app.misc.toBeEnabled = backup.misc.toBeEnabled;
                app.misc.disabled = backup.misc.disabled;
                app.misc.skip = backup.misc.skip;
                enableable = app.getEnableable;
                for e = 1:numel(enableable)
                    enableable(e).Enable = backup.enableable(e).Enable;
                end
                restoreText(app.load.currentTask, backup.load.currentTask);
                restoreText(app.load2.currentTask, backup.load2.currentTask);
                if ~app.misc.axToolbars
                    set(app.misc.menus, 'Enable', backup.misc.menus.Enable);
                end
            end
            
            drawnow
            
            function restoreButton(b, copyB)
                b.Visible = copyB.Visible;
                b.Enable = copyB.Enable;
            end
            
            function restoreText(t, copyT)
                t.Visible = copyT.Visible;
                if isprop(t, 'Enable'), t.Enable = copyT.Enable; end
                t.String = copyT.String;
                t.Position = copyT.Position;
            end
            
            function restoreAxes(a, copyA)
                a.Visible = copyA.Visible;
                a.XColor = copyA.XColor;
                a.YColor = copyA.YColor;
                a.ZColor = copyA.ZColor;
                a.Position = copyA.Position;
                a.View = copyA.View;
                a.XLimMode = 'auto';
                a.YLimMode = 'auto';
                a.ZLimMode = 'auto';
                if app.misc.axToolbars
                    a.Toolbar.Visible = copyA.Toolbar.Visible;
                end
            end
            
            function restoreAxesLimits(a, copyA)
                if strcmp(a.Visible, 'on')
                    axis(a, copyA.Limits);
                end
            end
            
            function p = restorePlot(p, copyP, varargin)
                if isempty(copyP), delete(p); p = gobjects(0); return, end
                if numel(p) < numel(copyP), delete(p); p = gobjects(numel(copyP),1); end
                if numel(p) > numel(copyP), delete(p(numel(copyP)+1:end)); p = p(1:numel(copyP)); end
                for i = 1:numel(copyP)
                    if ~ishghandle(p(i))
                        p(i) = plot(copyP(i).Parent, copyP(i).XData, copyP(i).YData,...
                            'Visible', copyP(i).Visible,...        
                            'HitTest', 'off',...
                            varargin{:});
                    else
                        p(i).Visible = copyP(i).Visible;
                        p(i).XData = copyP(i).XData;
                        p(i).YData = copyP(i).YData;
                        if isfield(copyP(i), 'ZData'), p(i).ZData = copyP(i).ZData; end
                    end
                end
            end
        end
        
        function fullReset2(app)
            % Reset tab 2 properties to exactly as they were after initialize().
            
            app.restoreBackupTab2(app.misc.initBackupTab2, 'fullReset');
        end
        
        function chooseXlsxCallback(app, ~, ~)
            % Called when user clicks load2.inputBtn. Allows user to pick a .xlsx file and
            % then updates load and analyze panels of tab 2.
            
            backup = app.storeBackupTab2;
            try
                % Let user pick .xlsx file.
                [name, path] = uigetfile('.xlsx');
                if name == 0
                    return
                end
                app.fullReset2;
                store = app.disable;
                drawnow;

                [~, ~, ext] = fileparts([path name]);
                if ~strcmp(ext, '.xlsx')
                    uiwait(errordlg('Input must be a .xlsx file.', 'Error Dialog', 'modal'))
                    app.restoreBackupTab2(backup);
                    return
                end

                % Set up variables.
                app.updateCurrentTask('Preparing to read data.');
                nSheets = app.numSheets([path name]);
                errorMsgs = cell(nSheets, 1);
                failed = false(nSheets, 1);
                pathsToTxts = cell(nSheets, 1);

                % Determine non-empty sheets in this .xlsx file.
                for i = 1:nSheets                
                    % For each sheet, save the data in cell A1 (the .txt path).
                    try
                        T = readtable([path name], 'Sheet', i, 'Range', 'A1:A1', 'ReadVariableNames', 0);
                        [txtpath, txtname] = fileparts(char(T{:,:}));
                        if isspace(txtname)
                            failed(i) = true; % Cell A1 is empty and does not reference a .txt file.
                            errorMsgs{i} = sprintf(['Sheet %d: Cell A1 should contain '...
                                'the path to a .txt file. This sheet has therefore been skipped.'], i);
                            continue
                        end
                        txtfile = fullfile(txtpath, [txtname '.txt']);
                        fileID = fopen(txtfile);
                        if fileID ~= -1
                            pathsToTxts{i} = txtfile;
                            fclose(fileID);
                        else
                            failed(i) = true; % There is a file referenced, but it cannot be opened - it may not exist.
                            errorMsgs{i} = sprintf(['Sheet %d: Unable to open "%s", the file '...
                                'referenced in cell A1. This sheet has therefore been skipped.'], i, txtfile);
                        end
                    catch ME
                        try
                            fclose(fileID);
                        catch
                        end
                        switch ME.identifier
                            case 'MATLAB:invalidConversion'
                                failed(i) = true; % There is no valid .txt file referenced in this sheet.
                                errorMsgs{i} = sprintf(['Sheet %d: Cell A1 should contain '...
                                    'the path to a .txt file. This sheet has therefore been skipped.'], i);
                            otherwise
                                rethrow(ME)
                        end
                    end
                end

                % Determine sampling frequency, events, event durations, step sizes, forces,
                % and bead positions.
                app.data2.events = cell(nSheets, 1);
                app.data2.durs = cell(nSheets, 1);
                app.data2.beadAStep1 = cell(nSheets, 1);
                app.data2.beadBStep1 = cell(nSheets, 1);
                app.data2.beadAStep2 = cell(nSheets, 1);
                app.data2.beadBStep2 = cell(nSheets, 1);
                app.data2.beadATotalStep = cell(nSheets, 1);
                app.data2.beadBTotalStep = cell(nSheets, 1);
                app.data2.forceA = cell(nSheets, 1);
                app.data2.forceB = cell(nSheets, 1);
                app.data2.beadA = cell(nSheets, 1);
                app.data2.beadB = cell(nSheets, 1);
                for i = 1:nSheets
                    app.updateCurrentTask(sprintf('Reading sheet %d.', i));

                    if failed(i)
                        % If this sheet does not reference a valid .txt file, ignore it.
                        continue
                    end

                    try
                        % Try reading data from xlsx file.
                        testNumeric = @(a, s) assert(isa(a, 'numeric') && ~isnan(a),...
                            'The %s was not a number.', s);
                        testOneMatch = @(a, s) assert(nnz(a)==1, ['There should be exactly '...
                            'one column containing the phrase ''%s'''], s);
                        
                        if verLessThan('matlab', '9.7')
                            firstCol = readtable([path name], 'Sheet', i, 'Range', 'A:A', 'ReadVariableNames', 0);
                            headerRows = find(contains(firstCol{:,:}, 'start (# pts)'))-1;
                            secondCol = table2cell(readtable([path name], 'Sheet', i, 'Range', ['B1:B' num2str(headerRows)], 'ReadVariableNames', 0));
                        else
                            opts = spreadsheetImportOptions('DataRange', 'A:A');
                            firstCol = readtable([path name], opts, 'Sheet', i, 'ReadVariableNames', 0);
                            headerRows = find(contains(firstCol{:,:}, 'start (# pts)'))-1;
                            opts = spreadsheetImportOptions('DataRange', ['B1:B' num2str(headerRows)]);
                            secondCol = table2cell(readtable([path name], opts, 'Sheet', i, 'ReadVariableNames', 0));
                        end
                        
                        % For each sheet, read the sampling frequency, number of events, and
                        % whether the beads should be flipped.
                        fs = secondCol{contains(firstCol{:,:}, 'sampling frequency (Hz)')};
                        if ischar(fs)
                            fs = str2double(fs);
                        end
                        testNumeric(fs, 'sampling frequency');
                        
                        N = secondCol{contains(firstCol{:,:}, 'number of events')};
                        if ischar(N)
                            N = str2double(N);
                        end
                        testNumeric(N, 'number of events');
                        
                        flip = secondCol{contains(firstCol{:,:}, 'flipped')};
                        if ischar(flip)
                            assert(strcmpi(flip,'true')||strcmpi(flip,'false')...
                                ||strcmpi(flip,'1')||strcmpi(flip,'0'),...
                                'The value for flipped was neither true nor false.');
                            flip = strcmpi(flip,'true') || strcmpi(flip,'1');
                        end
                        assert(flip==1||flip==0||flip==true||flip==false,...
                            'The value for flipped was neither true nor false.');

                        % Also read the table on this sheet.
                        colNames = table2cell(readtable([path name], 'Sheet', i, 'Range', [num2str(headerRows+1) ':' num2str(headerRows+1)], 'ReadVariableNames', 0));
                        colNames(~cellfun(@ischar,colNames)) = {''};
                        xlsxData = table2array(readtable([path name], 'Sheet', i, 'Range', [num2str(headerRows+2) ':' num2str(headerRows+1+N)], 'ReadVariableNames', 0));
                        testNumeric = @(a, s) assert(isa(a, 'numeric') && ~any(any(isnan(a))),...
                            'The %s contained non-numeric values.', s);
                        
                        startIdx = contains(colNames, 'start (# pts)');
                        stopIdx = contains(colNames, 'end (# pts)');
                        durIdx = contains(colNames, 'duration (# pts)');
                        Astep1Idx = contains(colNames, 'A step 1');
                        Bstep1Idx = contains(colNames, 'B step 1');
                        Astep2Idx = contains(colNames, 'A step 2');
                        Bstep2Idx = contains(colNames, 'B step 2');
                        AtotalStepIdx = contains(colNames, 'A total step');
                        BtotalStepIdx = contains(colNames, 'B total step');
                        AforceDuringIdx = contains(colNames, 'A avg force during');
                        AforceAfterIdx = contains(colNames, 'A force after detach');
                        BforceDuringIdx = contains(colNames, 'B avg force during');
                        BforceAfterIdx = contains(colNames, 'B force after detach');
                        
                        testOneMatch(startIdx, 'start (# pts)');
                        testOneMatch(stopIdx, 'end (# pts)');
                        testOneMatch(durIdx, 'duration (# pts)');
                        testOneMatch(Astep1Idx, 'A step 1');
                        testOneMatch(Bstep1Idx, 'B step 1');
                        testOneMatch(Astep2Idx, 'A step 2');
                        testOneMatch(Bstep2Idx, 'B step 2');
                        testOneMatch(AtotalStepIdx, 'A total step');
                        testOneMatch(BtotalStepIdx, 'B total step');
                        testOneMatch(AforceDuringIdx, 'A avg force during');
                        testOneMatch(AforceAfterIdx, 'A force after detach');
                        testOneMatch(BforceDuringIdx, 'B avg force during');
                        testOneMatch(BforceAfterIdx, 'B force after detach');
                        
                        testNumeric(xlsxData(:,startIdx | stopIdx), 'events');
                        testNumeric(xlsxData(:,durIdx), 'event durations');
                        testNumeric(xlsxData(:,Astep1Idx), 'bead A step 1 sizes');
                        testNumeric(xlsxData(:,Bstep1Idx), 'bead B step 1 sizes');
                        testNumeric(xlsxData(:,Astep2Idx), 'bead A step 2 sizes');
                        testNumeric(xlsxData(:,Bstep2Idx), 'bead B step 2 sizes');
                        testNumeric(xlsxData(:,AtotalStepIdx), 'bead A total step sizes');
                        testNumeric(xlsxData(:,BtotalStepIdx), 'bead B total step sizes');
                        testNumeric(xlsxData(:,AforceDuringIdx | AforceAfterIdx), 'bead A forces');
                        testNumeric(xlsxData(:,BforceDuringIdx | BforceAfterIdx), 'bead B forces');
                        
                        % If all values are valid, store them.
                        if isempty(app.data2.fs)
                            app.data2.fs = fs;
                        end
                        app.data2.events{i} = xlsxData(:,startIdx | stopIdx);
                        app.data2.durs{i} = xlsxData(:,durIdx);
                        app.data2.beadAStep1{i} = xlsxData(:,Astep1Idx);
                        app.data2.beadBStep1{i} = xlsxData(:,Bstep1Idx);
                        app.data2.beadAStep2{i} = xlsxData(:,Astep2Idx);
                        app.data2.beadBStep2{i} = xlsxData(:,Bstep2Idx);
                        app.data2.beadATotalStep{i} = xlsxData(:,AtotalStepIdx);
                        app.data2.beadBTotalStep{i} = xlsxData(:,BtotalStepIdx);
                        app.data2.forceA{i} = xlsxData(:,AforceDuringIdx)-xlsxData(:,AforceAfterIdx);
                        app.data2.forceB{i} = xlsxData(:,BforceDuringIdx)-xlsxData(:,BforceAfterIdx);
                    catch ME
                        % If any values are invalid, skip the entire sheet.
                        errorMsgs{i} = sprintf('Sheet %d: %s This sheet has therefore been skipped.', i, ME.message);
                        failed(i) = true;
                        continue
                    end

                    txtfile = pathsToTxts{i};
                    app.updateCurrentTask(sprintf('Reading %s.', txtfile));
                    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames') % Suppress warning output.
                    try
                        % Try reading data from txt file.
                        try
                            % Try using the data input function.
                            [~, ~, ~, CALA, CALB, txtData, ~] = app.load.txtInputFcn(txtfile);
                        catch ME
                            % If the file could not be read, continue to the next sheet.
                            errorMsgs{i} = sprintf(['Sheet %d: There was an error while reading "%s":\n\n'...
                                '%s\n\nThis sheet has therefore been skipped.'],...
                                i, txtfile, ME.message);
                            failed(i) = true;
                            continue
                        end
                        % If file was successfully read, now test the outputs.
                        testNumeric = @(a, s) assert(isa(a, 'numeric') && ~isnan(a),...
                            'chooseXlsxCallback:outputFailed',...
                            'the %s was not a number.', s);
                        testNumeric(CALA, 'value of CAL for bead A');
                        testNumeric(CALB, 'value of CAL for bead B');
                        assert(istable(txtData),...
                            'chooseXlsxCallback:outputFailed',...
                            'the final variable returned by the data input function had type %s but should be a table.',...
                            class(txtData));
                        j = ismember({'BeadAPos','BeadBPos'}, txtData.Properties.VariableNames);
                        assert(j(1),...
                            'chooseXlsxCallback:outputFailed',...
                            'the table returned by the data input function did not contain a column labeled "BeadAPos".');
                        assert(j(2),...
                            'chooseXlsxCallback:outputFailed',...
                            'the table returned by the data input function did not contain a column labeled "BeadBPos".');
                    catch ME
                        if strcmp(ME.identifier, 'chooseXlsxCallback:outputFailed')
                            % If any of the outputs are invalid:
                            errorMsgs{i} = sprintf(['Sheet %d: When trying to read "%s", %s '...
                                'This sheet has therefore been skipped.'],...
                                i, txtfile, ME.message);
                        else
                            % If there was some other error:
                            errorMsgs{i} = sprintf(['Sheet %d: There was an unexpected error when '...
                                'trying to read "%s":\n\n%s\n\nThis sheet has therefore been skipped.'],...
                                i, txtfile, ME.message);
                        end
                        failed(i) = true;
                        continue
                    end
                    warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')

                    % Flip if needed.
                    if flip
                        CALA = -CALA;
                        CALB = -CALB;
                    end

                    % Convert from V to nm.
                    BeadAPos = txtData.BeadAPos * CALA;
                    BeadBPos = txtData.BeadBPos * CALB;

                    try
                        % Determine which data should be removed.
                        deselectIndices = true(length(BeadAPos), 1); % Used to determine which indices of averageBeads user opted to remove via deselect.
                        [col, j] = app.getExcelCol(2); % [second excel column as a string (so, 'B'), index of next excel column (here, 3)]
                        row = find(contains(firstCol{:,:}, 'ignore data between these indices (start 1, stop 1, start 2, stop 2, etc)'));
                        res = table2array(readtable([path name], 'Sheet', i, 'Range', [col num2str(row) ':' col num2str(row)], 'ReadVariableNames', 0));
                        while ~isnan(res)
                            % As long as the last cell read has data, update deselectIndices and
                            % read the next cell in row 5.
                            assert(isnumeric(res) && ~isnan(res),...
                                'chooseXlsxCallback:badIndices',...
                                'at least one of the indices marking data to ignore was not a number.');
                            assert(res >= 1 && res <= length(BeadAPos),...
                                'chooseXlsxCallback:badIndices',...
                                ['at least one of the indices marking data to ignore was out of bounds. Each '...
                                'index needs to be between 1 and the length of the data, %d.'],...
                                length(BeadAPos));
                            deselectIndices(res) = false;
                            [col, j] = app.getExcelCol(j); % [current excel column as a string, index of next excel column]
                            res = table2array(readtable([path name], 'Sheet', i, 'Range', [col num2str(row) ':' col num2str(row)], 'ReadVariableNames', 0));
                        end
                        ends = find(~deselectIndices); % Edges of patches of data to be removed.
                        for j = 1:2:length(ends)-1
                            % For each patch, set data between the edges to be removed.
                            deselectIndices(ends(j):ends(j+1)) = false;
                        end
                        % Remove that data.
                        app.data2.beadA{i} = BeadAPos(deselectIndices);
                        app.data2.beadB{i} = BeadBPos(deselectIndices);
                    catch ME
                        if strcmp(ME.identifier, 'chooseXlsxCallback:badIndices')
                            % If there was an error reading the deselect indices:
                            errorMsgs{i} = sprintf(['Sheet %d: When trying to read "%s", %s '...
                                'This sheet has therefore been skipped.'], i, txtfile, ME.message);
                        else
                            % If there was some other error:
                            errorMsgs{i} = sprintf(['Sheet %d: There was an unexpected error when '...
                                'trying to read "%s":\n\n%s\n\nThis sheet has therefore been skipped.'],...
                                i, txtfile, ME.message);
                        end
                        failed(i) = true;
                        continue
                    end
                end

                % Inform user if any .xlsx files failed to load properly.
                app.updateCurrentTask('Cleaning up.');
                if any(failed)
                    if all(failed)
                        uiwait(warndlg(sprintf(['None of the sheets in your .xlsx file could be read. '...
                            'Below is the error message for the first sheet:\n\n\n%s'], errorMsgs{1}),...
                            'Error Dialog', 'modal'))
                        app.restoreBackupTab2(backup);
                        return
                    else
                        app.data2.events = app.data2.events(~failed);
                        app.data2.durs = app.data2.durs(~failed);
                        app.data2.beadAStep1 = app.data2.beadAStep1(~failed);
                        app.data2.beadBStep1 = app.data2.beadBStep1(~failed);
                        app.data2.beadAStep2 = app.data2.beadAStep2(~failed);
                        app.data2.beadBStep2 = app.data2.beadBStep2(~failed);
                        app.data2.beadATotalStep = app.data2.beadATotalStep(~failed);
                        app.data2.beadBTotalStep = app.data2.beadBTotalStep(~failed);
                        app.data2.forceA = app.data2.forceA(~failed);
                        app.data2.forceB = app.data2.forceB(~failed);
                        app.data2.beadA = app.data2.beadA(~failed);
                        app.data2.beadB = app.data2.beadB(~failed);

                        warndlg(sprintf(['The program encountered the following error(s) when '...
                            'reading your .xlsx file:\n\n\n%s'],...
                            strip(sprintf('%s\n\n\n', errorMsgs{failed}))))
                    end
                end

                % Update the table.
                app.updateCurrentTask('Filling the table.');
                N = cellfun(@(x) size(x, 1), app.data2.events);
                app.load2.tableData = ...
                    [cellstr(num2str(find(~failed))),... % The sheet number.
                    pathsToTxts(~failed),...             % The .txt filename.
                    cellstr(num2str(N)),...              % The number of events.
                    num2cell(true(numel(N),1))];         % Whether to include this file in the analysis.
                app.load2.tWidths = app.fillTable(app.load2.table, app.load2.tableData, app.load2.maxHeight);
                app.load2.table.CellEditCallback = @app.cellEditCallback;
                drawnow;
                
                % Plot misc axes.
                app.plotMisc2;
                
                % Calculate and plot ensemble averages.
                [app.analyze2.ensemble.start, app.analyze2.ensemble.end] = app.calcEnsemble2;
                app.plotEnsemble2;
                if app.misc.optim
                    [app.analyze2.ensemble.startFit, app.analyze2.ensemble.endFit,...
                        app.analyze2.ensemble.p] = app.calcEnsembleFits2;
                    app.plotEnsembleFits2;
                end
                
                app.load2.whichBeadBtnGrp.Visible = 'on';
                app.load2.beadABtn.Visible = 'on';
                app.load2.beadBBtn.Visible = 'on';
                app.load2.bothBeadsBtn.Visible = 'on';
                
                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab2(backup);
            end
        end
        
        function steps = getStepSizes2(app, bead)
            % For each event, for the requested bead, returns the step 1, step 2, and
            % total step sizes, in a cell array.
            
            if strcmp(bead, 'B')
                steps = [app.data2.beadBStep1 app.data2.beadBStep2 app.data2.beadBTotalStep];
            else
                steps = [app.data2.beadAStep1 app.data2.beadAStep2 app.data2.beadATotalStep];
            end
        end
        
        function cellEditCallback(app, ~, eventData)
            % Called when user changes the Include value of any file in the table in tab
            % 2.
            
            backup = app.storeBackupTab2;
            try
                store = app.disable;
                drawnow;

                use = cellfun(@(c) c, app.load2.table.Data(:,end));
                if all(~use)
                    % If the user tries to uncheck the last row, do not let them.
                    app.load2.table.Data{eventData.Indices(1),end} = true;
                else
                    % Otherwise, update the plots.
                    app.plotMisc2;
                    [app.analyze2.ensemble.start, app.analyze2.ensemble.end] = app.calcEnsemble2;
                    app.plotEnsemble2;
                    if app.misc.optim
                        [app.analyze2.ensemble.startFit, app.analyze2.ensemble.endFit,...
                            app.analyze2.ensemble.p] = app.calcEnsembleFits2;
                        app.plotEnsembleFits2;
                    end
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab2(backup);
            end
        end
        
        function whichBeadCallback2(app, src, eventData)
            % Called when user changes which bead(s) should be analyzed.
            
            backup = app.storeBackupTab2;
            try
                store = app.disable;
                drawnow;

                if any(strcmp(app.analyze2.misc.data, {'Step 1 Size', 'Step 2 Size', 'Total Step Size',...
                        'Event Duration / Force on Selected Bead'}))
                    % Step size plots are dependent on the bead.
                    app.plotMisc2;
                end

                [app.analyze2.ensemble.start, app.analyze2.ensemble.end] = app.calcEnsemble2;
                app.plotEnsemble2;
                if app.misc.optim
                    [app.analyze2.ensemble.startFit, app.analyze2.ensemble.endFit,...
                        app.analyze2.ensemble.p] = app.calcEnsembleFits2;
                    app.plotEnsembleFits2;
                end

                app.restore(store);
            catch ME
                if ~isvalid(app)
                    return
                end
                app.unexpectedError(ME);
                app.restoreBackupTab2(backup);
                src.SelectedObject = eventData.OldValue;
            end
        end
        
        function plotMisc2(app)
            % Calculates and plots various data and fits on the misc axes in analyze panel
            % of tab 2. Uses value of app.analyze2.misc.data to determine data and, if
            % appropriate, fit. Defaults to event durations if needed.
            
            use = cellfun(@(c) c, app.load2.table.Data(:,end));
            
            if strcmp(app.analyze2.misc.data, 'Event Duration / Force on Selected Bead') &&...
                    app.load2.bothBeadsBtn.Value
                app.analyze2.misc.data = 'Event Duration';
            end
            
            % Determine which data and distribution to use.
            ylab = 'cumulative distribution';
            switch app.analyze2.misc.data
                case 'Step 1 Size'
                    stepA = app.getStepSizes2('A');
                    stepB = app.getStepSizes2('B');
                    if app.load2.beadABtn.Value
                        % Just bead A.
                        Y = cell2mat(stepA(use,1));
                        xlab = 'step 1 size (nm), estimated by bead A';
                    elseif app.load2.beadBBtn.Value
                        % Just bead B.
                        Y = cell2mat(stepB(use,1));
                        xlab = 'step 1 size (nm), estimated by bead B';
                    else
                        % Both.
                        Y = (cell2mat(stepA(use,1))+cell2mat(stepB(use,1)))/2;
                        xlab = 'step 1 size (nm), estimated by avg of both beads';
                    end
                    type = 'normcdf';
                    app.updateCurrentTask('Plotting the distribution of substep 1 sizes.');
                case 'Step 2 Size'
                    stepA = app.getStepSizes2('A');
                    stepB = app.getStepSizes2('B');
                    if app.load2.beadABtn.Value
                        % Just bead A.
                        Y = cell2mat(stepA(use,2));
                        xlab = 'step 2 size (nm), estimated by bead A';
                    elseif app.load2.beadBBtn.Value
                        % Just bead B.
                        Y = cell2mat(stepB(use,2));
                        xlab = 'step 2 size (nm), estimated by bead B';
                    else
                        % Both.
                        Y = (cell2mat(stepA(use,2))+cell2mat(stepB(use,2)))/2;
                        xlab = 'step 2 size (nm), estimated by avg of both beads';
                    end
                    type = 'normcdf';
                    app.updateCurrentTask('Plotting the distribution of substep 2 sizes.');
                case 'Total Step Size'
                    stepA = app.getStepSizes2('A');
                    stepB = app.getStepSizes2('B');
                    if app.load2.beadABtn.Value
                        % Just bead A.
                        Y = cell2mat(stepA(use,3));
                        xlab = 'total step size (nm), estimated by bead A';
                    elseif app.load2.beadBBtn.Value
                        % Just bead B.
                        Y = cell2mat(stepB(use,3));
                        xlab = 'total step size (nm), estimated by bead B';
                    else
                        % Both.
                        Y = (cell2mat(stepA(use,3))+cell2mat(stepB(use,3)))/2;
                        xlab = 'total step size (nm), estimated by avg of both beads';
                    end
                    type = 'normcdf';
                    app.updateCurrentTask('Plotting the distribution of total step sizes.');
                case 'Event Duration / Force on Selected Bead'
                    if app.load2.beadABtn.Value
                        bead = 'A';
                        X = cell2mat(app.data2.forceA(use));
                    else
                        bead = 'B';
                        X = cell2mat(app.data2.forceB(use));
                    end
                    Y = cell2mat(app.data2.durs(use))/app.data2.fs;
                    type = 'semilogy';
                    ylab = 'event duration (s)';
                    xlab = ['force on bead ' bead ' (pN)'];
                    app.updateCurrentTask(['Plotting the event durations over the force on bead ' bead '.']);
                otherwise
                    Y = cell2mat(app.data2.durs(use))/app.data2.fs;
                    type = 'expcdf';
                    xlab = 'event duration (s)';
                    app.updateCurrentTask('Plotting the distribution of event durations.');
            end
            
            % Determine the data and, if appropriate, the fit to plot.
            switch type
                case 'expcdf'
                    pdf = @(x, k) k*exp(-k*x); % Exponential distribution probability density function.
                    mle = @(k) -sum(log(pdf(Y, k))); % Function to minimize during maximum likelihood estimation.
                    p = fminsearch(mle, 0.0005); % k value which minimizes the mle function.
                    fit = @(x) 1 - exp(-p*x); % Fit to exponential cumulative distribution.
                    n = fit(min(Y))*(length(Y)); % Estimated number of events that are shorter than smallest measured event.
                    X = linspace(min(Y), max(Y), 1000); % Range of D for the cumulative distribution.
                    fitX = linspace(0, max(Y), 1000); % Range of D for the fit.
                    Y = (sum(Y<=X,1)+n) / (length(Y)+n); % Cumulative distribution function.
                    fitY = fit(fitX);
                    s = ['k = ' sprintf('%.3g', p) ' s^{-1}'];
                    pos = [0.5 0.5];
                case 'normcdf'
                    pdf = @(x, p) (1/(p(2)*sqrt(2*pi)))*exp(-0.5*((x-p(1))/p(2)).^2); % Normal distribution probability density function.
                    mle = @(p) -sum(log(pdf(Y, p))); % Function to minimize during maximum likelihood estimation.
                    p = fminsearch(mle, [mean(Y) std(Y)]); % p value which minimizes the mle function.
                    fit = @(x) 0.5*erfc((p(1)-x)/(p(2)*sqrt(2))); % Fit to normal cumulative distribution.
                    X = linspace(min(Y), max(Y), 1000); % Range of D for the cumulative distribution.
                    fitX = p(1)+p(2)*-sqrt(2)*erfcinv(2*linspace(0.001,0.999,1000)); % Range of D for the fit.
                    Y = sum(Y<=X,1) / length(Y); % Cumulative distribution function.
                    fitY = fit(fitX);
                    s = sprintf('\\mu = %.3g nm\\newline\\sigma = %.3g nm', [p(1) p(2)]);
                    pos = [0.25, 5];
                case 'semilogy'
                    fitX = [];
                    fitY = [];
                    s = '';
                    pos = [0, 0];
            end
            
            % Update plots.
            app.changeData(app.analyze2.misc.miscPlot, X, Y);
            app.changeData(app.analyze2.misc.miscFitPlot, fitX, fitY);
            xlabel(app.analyze2.misc.miscAxes, xlab);
            ylabel(app.analyze2.misc.miscAxes, ylab);
            
            % If scatter plot, change marker, line style, y scale, axis limits.
            if strcmp(type, 'semilogy')
                app.analyze2.misc.miscPlot.Marker = 'o';
                app.analyze2.misc.miscPlot.LineStyle = 'none';
                app.analyze2.misc.miscAxes.YScale = 'log';
                axis(app.analyze2.misc.miscAxes, 'tight');
                limx = app.analyze2.misc.miscAxes.XLim;
                limy = app.analyze2.misc.miscAxes.YLim;
                app.analyze2.misc.miscAxes.XLim = limx + [-0.1 0.1]*(limx(2)-limx(1));
                app.analyze2.misc.miscAxes.YLim = 10.^(log10(limy) + [-0.1 0.1]*(log10(limy(2))-log10(limy(1))));
            else
                app.analyze2.misc.miscPlot.Marker = 'none';
                app.analyze2.misc.miscPlot.LineStyle = '-';
                app.analyze2.misc.miscAxes.YScale = 'linear';
                axis(app.analyze2.misc.miscAxes, 'tight');
            end
            
            % Update label.
            app.analyze2.misc.miscLabel.String = '';
            [~, x, y] = app.placeText(app.analyze2.misc.miscAxes, s, pos, app.misc.fontSize);
            app.analyze2.misc.miscLabel.Position = [x y];
            app.analyze2.misc.miscLabel.String = s;
            app.analyze2.misc.miscLabel.Visible = 'on';
            
            % Show axes and toolbar.
            app.analyze2.misc.miscAxes.Visible = 'on';
            
            % Resize axes based on yticklabels and ylabel.
            ti = app.analyze2.misc.miscAxes.TightInset;
            pos = app.analyze2.misc.miscAxes.Position;
            pos(3) = 0.495 - pos(1) - ti(3);
            app.analyze2.misc.miscAxes.Position = pos;
            
            % Update stored parameters.
            switch app.analyze2.misc.data
                case 'Step 1 Size'
                    app.analyze2.misc.meanStep1 = p(1);
                    app.analyze2.misc.stdStep1 = p(2);
                case 'Step 2 Size'
                    app.analyze2.misc.meanStep2 = p(1);
                    app.analyze2.misc.stdStep2 = p(2);
                case 'Total Step Size'
                    app.analyze2.misc.meanTotalStep = p(1);
                    app.analyze2.misc.stdTotalStep = p(2);
                case 'Event Duration / Force on Selected Bead'
                otherwise
                    app.analyze2.misc.kDur = p;
            end
            
            drawnow;
        end
        
        function [avgbeg, avgend] = calcEnsemble2(app)
            % Calculate and plot the ensemble averages and fits on the ensemble axes in
            % analyze panel of tab2.
            
            app.updateCurrentTask('Calculating the ensemble averages.');
            use = cellfun(@(c) c, app.load2.table.Data(:,end));
            
            durs = app.data2.durs(use);
            A = app.data2.beadA(use);
            B = app.data2.beadB(use);
            events = app.data2.events(use);
            
            if app.load2.beadABtn.Value
                % Just bead A.
                averageBeads = A;
            elseif app.load2.beadBBtn.Value
                % Just bead B.
                averageBeads = B;
            else
                % Both.
                averageBeads = cell(size(A));
                for i = 1:numel(A)
                    averageBeads{i} = (A{i}+B{i})/2;
                end
            end
            
            pts_bef = round(app.analyze.ensemble.pts_bef * app.data2.fs);
            pts_after = round(app.analyze.ensemble.pts_after * app.data2.fs);
            exten_pos = round(app.analyze.ensemble.exten_pos * app.data2.fs);
            
            % Initialize arrays to hold event beginnings and endings.
            maxDur = max(cell2mat(durs));
            N = sum(cellfun(@(x) size(x, 1), events));
            beginnings = ones(maxDur + pts_bef - exten_pos, N);
            endings = ones(maxDur + pts_after - exten_pos, N);
            
            k = 0;
            for i = 1:length(averageBeads)
                % For each file:
                
                % Get the beads, events, and durations for that file.
                fileAverageBeads = averageBeads{i};
                fileDurs = durs{i};
                fileEvents = events{i};
                N = size(fileEvents, 1);
                
                % Calculate average of event beginnings.
                beginning = ones(maxDur + pts_bef - exten_pos, size(fileEvents, 1));
                for j = 1:size(fileEvents, 1)
                    % For each event...
                    if fileEvents(j,1) - pts_bef < 1
                        % ...if event window starts N points before start of file, then
                        % first N values of event window will match first value in file...
                        event = [fileAverageBeads(1)*ones(pts_bef - fileEvents(j,1),1); ...
                            fileAverageBeads(1:fileEvents(j,2) - exten_pos)];
                    else
                        % ...otherwise determine the event window...
                        event = fileAverageBeads(fileEvents(j,1) + 1 - pts_bef:fileEvents(j,2) - exten_pos);
                    end
                    % ...and calculate the extension so that each event beginning is the
                    % same size as other event beginnings.
                    extension = mean(fileAverageBeads(fileEvents(j,2) - exten_pos - app.data2.fs/200:fileEvents(j,2) - exten_pos))*ones(maxDur - fileDurs(j),1);
                    beginning(:,j) = [event; extension];
                end
                beginnings(:, k+1:k+N) = beginning; % Store the event beginnings.
                
                % Calculate average of event endings.
                ending = ones(maxDur + pts_after - exten_pos, size(fileEvents, 1));
                for j = 1:size(fileEvents, 1)
                    % For each event...
                    if fileEvents(j,2) + pts_after > length(fileAverageBeads)
                        % ...if event window ends N points after end of file, then last N
                        % values of event window will match last value in file...
                        event = [fileAverageBeads(fileEvents(j,1) + 1 + exten_pos:end); ...
                            fileAverageBeads(end)*ones(pts_after + fileEvents(j,2) - length(fileAverageBeads),1)];
                    else
                        % ...otherwise determine the event window...
                        event = fileAverageBeads(fileEvents(j,1) + 1 + exten_pos:fileEvents(j,2) + pts_after);
                    end
                    % ...and calculate the extension so that each event ending is the same
                    % size as other event endings.
                    extension = mean(fileAverageBeads(fileEvents(j,1) + 1 + exten_pos:fileEvents(j,1) + 1 + exten_pos + app.data2.fs/200))*ones(maxDur - fileDurs(j),1);
                    ending(:,j) = [extension; event];
                end            
                endings(:, k+1:k+N) = ending; % Store the event endings.
                k = k+N;
            end
            
            % Average the event beginnings and endings.
            avgbeg = mean(beginnings, 2);
            avgend = mean(endings, 2);
            
            % Adjust averages so beginning average starts at 0.
            offset = movmean(avgbeg(1:pts_bef), app.data2.fs/10);
            avgbeg = avgbeg - offset(end);
            avgend = avgend - offset(end);
        end
        
        function plotEnsemble2(app)
            % Plots the combined ensemble averages on the ensemble axes in analyze panel
            % of tab 2.
            
            app.updateCurrentTask('Plotting the ensemble averages.');
            
            startAvg = app.analyze2.ensemble.start;
            endAvg = app.analyze2.ensemble.end;
            
            if app.load2.beadABtn.Value
            	ylab = 'position of bead A (nm)';
            elseif app.load2.beadBBtn.Value
                ylab = 'position of bead B (nm)';
            else
                ylab = 'avg pos of both beads (nm)';
            end
            
            % Calculate time ranges for each average.
            startTime = (1:length(startAvg)) / app.data2.fs;
            endTime = (length(startTime)+(1:length(endAvg))) / app.data2.fs;
            
            % Update plots.
            app.changeData(app.analyze2.ensemble.startPlot, startTime, startAvg);
            app.changeData(app.analyze2.ensemble.endPlot, endTime, endAvg);
            axis(app.analyze2.ensemble.ensembleAxes, 'tight');
            app.analyze2.ensemble.ensembleAxes.Visible = 'on';
            ylabel(app.analyze2.ensemble.ensembleAxes, ylab);
            
            % Resize axes based on yticklabels and ylabel.
            ti = app.analyze2.ensemble.ensembleAxes.TightInset;
            pos = app.analyze2.ensemble.ensembleAxes.Position;
            pos(3) = 0.99 - pos(1) - ti(3);
            app.analyze2.ensemble.ensembleAxes.Position = pos;
            
            % Hide fits and labels.
            app.analyze2.ensemble.startFitPlot.Visible = 'off';
            app.analyze2.ensemble.endFitPlot.Visible = 'off';
            app.analyze2.ensemble.startLabel.Visible = 'off';
            app.analyze2.ensemble.endLabel.Visible = 'off';
            app.analyze2.ensemble.stepLabel.Visible = 'off';
            
            drawnow;
        end
        
        function [startFit, endFit, p] = calcEnsembleFits2(app)
            % Calculates fits of the combined ensemble averages.
            
            app.updateCurrentTask('Fitting the ensemble averages.');
            
            fs = app.data2.fs;
            pts_bef = round(app.analyze.ensemble.pts_bef * fs);
            pts_after = round(app.analyze.ensemble.pts_after * fs);
            skip = round(app.analyze.ensemble.num_pts_to_skip);
            
            % Isolate the exponential portions of each average.
            avgbegexp = app.analyze2.ensemble.start(pts_bef+1+skip:end);
            startTime = (skip-1+(1:length(avgbegexp))) / fs;
            avgendexp = app.analyze2.ensemble.end(1:end-pts_after-skip);
            endTime = ((-length(avgendexp):-1)+1-skip) / fs;
            
            % Initialize variables.
            p = [];
            startFit = NaN(1, length(avgbegexp));
            endFit = NaN(1, length(avgendexp));
            t = timer('StartFcn', @(varargin) set(app.load2.skipBtn, 'Enable', 'on'),...
                'TimerFcn', @(varargin) set(app.load2.skipBtn, 'Visible', 'on'),...
                'StartDelay', 3);
            
            app.updateCurrentTask('Fitting the ensemble averages.');
            try
                start(t)
                [startFit, endFit, p] = app.fitexp(avgbegexp, startTime, avgendexp, endTime, app.analyze2.ensemble.globalFit); % [start fit, end fit, fit parameters]
                stop(t)
                app.load2.skipBtn.Visible = 'off';
            catch
                stop(t)
                app.load2.skipBtn.Visible = 'off';
            end
        end
        
        function plotEnsembleFits2(app)
            % Plots the combined ensemble average fits and labels on the ensemble axes in
            % analyze panel of tab 2.
            
            app.updateCurrentTask('Plotting the ensemble average fits.');
            
            startFit = app.analyze2.ensemble.startFit;
            endFit = app.analyze2.ensemble.endFit;
            pts_bef = round(app.analyze.ensemble.pts_bef * app.data2.fs);
            skip = round(app.analyze.ensemble.num_pts_to_skip);
            
            app.analyze2.ensemble.startLabel.String = '';
            app.analyze2.ensemble.endLabel.String = '';
            app.analyze2.ensemble.stepLabel.String = '';
            
            if ~any(isnan(startFit)) && ~isempty(app.analyze2.ensemble.p)
                % If start fit exists, update the start fit plot...
                startTime = (pts_bef+1+skip+(1:length(startFit))) / app.data2.fs;
                app.changeData(app.analyze2.ensemble.startFitPlot, startTime, startFit);
                
                % ...and label.
                kf = app.analyze2.ensemble.p(1);
                s = sprintf('k_f = %.3g s^{-1}', kf);
                [~, x, y] = app.placeText(app.analyze2.ensemble.ensembleAxes, s, [0.25, 0.5], app.misc.fontSize);
                app.analyze2.ensemble.startLabel.Position = [x y];
                app.analyze2.ensemble.startLabel.String = s;
                app.analyze2.ensemble.startLabel.Visible = 'on';
            end
            
            if ~any(isnan(endFit)) && ~isempty(app.analyze2.ensemble.p)
                % If end fit exists, update the end fit plot...
                endTime = (pts_bef+1+skip+length(startFit)+(1:length(endFit))) / app.data2.fs;
                app.changeData(app.analyze2.ensemble.endFitPlot, endTime, endFit);
                
                % ...and label.
                kr = app.analyze2.ensemble.p(2);
                s = sprintf('k_r = %.3g s^{-1}', kr);
                [~, x, y] = app.placeText(app.analyze2.ensemble.ensembleAxes, s, [0.75, 0.5], app.misc.fontSize);
                app.analyze2.ensemble.endLabel.Position = [x y];
                app.analyze2.ensemble.endLabel.String = s;
                app.analyze2.ensemble.endLabel.Visible = 'on';
            end
            
            if ~isempty(app.analyze2.ensemble.p)
                % Update step labels.
                step1 = app.analyze2.ensemble.p(3);
                totalStep = app.analyze2.ensemble.p(4);
                s = sprintf('step 1: %.3g nm\ntotal step: %.3g nm', step1, totalStep);
                limy = ylim(app.analyze2.ensemble.ensembleAxes);
                [~, x, y] = app.placeText(app.analyze2.ensemble.ensembleAxes, s, [0.5, -limy(1)/(limy(2)-limy(1))], app.misc.fontSize);
                app.analyze2.ensemble.stepLabel.Position = [x y];
                app.analyze2.ensemble.stepLabel.String = s;
                app.analyze2.ensemble.stepLabel.Visible = 'on';
            end
            
            axis(app.analyze2.ensemble.ensembleAxes, 'tight');
            
            drawnow;
        end
    end
    
    % Private static methods
    methods (Static, Access = private)
        
        function sHeader = defaultTxtOutput(fs, KA, KB, CALA, CALB, tData, header)
            % DEFAULTTXTOUTPUT Formats the header and column names for a text file.
            %   sHeader = DEFAULTTXTOUTPUT(fs, KA, KB, CALA, CALB, tData, header) accepts
            %   the variables fs, KA, KB, CALA, CALB, tData, and header exactly as they
            %   were returned by the app's .txt input function.
            %
            %   The output sHeader is a string containing the header and column names.
            %   sHeader will be printed directly to the text file, and so sHeader is
            %   formatted accordingly, with tabs (\t) separating columns of text and
            %   newline characters (\n) separating rows of text. The columns for beads A
            %   and B are renamed from 'BeadAPos' and 'BeadBPos' to 'Trap1X' and 'Trap2X'.
            
            try
                % First assume header is formatted as the output from defaultTxtInput()
                % (the app's default txtInputFcn).
            	assert(isnumeric(str2double(header{2}{strcmp(header{1}, 'K1')})));
                assert(isnumeric(str2double(header{2}{strcmp(header{1}, 'K3')})));
                assert(isnumeric(str2double(header{2}{strcmp(header{1}, 'CAL1')})));
                assert(isnumeric(str2double(header{2}{strcmp(header{1}, 'CAL3')})));
                assert(isnumeric(str2double(header{2}{strcmp(header{1}, 'Sample Rate')})));
                prepHeader = [header{1}(:)'; header{2}(:)'];
            catch
                % If header is not formatted in that way, then user must have supplied a
                % txtInputFcn but not a txtOutputFcn. Use the provided fs, KA,
                % KB, CALA, and CALB to create a header in the default format.
                header = cell(2,1);
                header{1} = {'Sample Rate', 'K1', 'K3', 'CAL1', 'CAL3'};
                header{2} = split(num2str([fs KA KB CALA CALB]));
                prepHeader = [header{1}(:)'; header{2}(:)'];
            end
            
            % Rename the columns that contain the data.
            tData.Properties.VariableNames{'BeadAPos'} = 'Trap1X';
            tData.Properties.VariableNames{'BeadBPos'} = 'Trap2X';
            
            % Format sHeader using the supplied header and the column names from tData.
            colNames = tData.Properties.VariableNames;
            sHeader = [sprintf('%s\t%s\n', prepHeader{:})...
                sprintf('%s\t', colNames{1:end-1})...
                sprintf('%s\n', colNames{end})];
        end
        
        function [norm, xco, yco] = placeText(ax, s, target, fs)
            % PLACETEXT Determines optimal placement of text on a set of axes.
            %   norm = PLACETEXT(ax, string, target, fs) returns the normalized X and Y
            %   coordinates for a hypothetical Text object t with String s and FontSize fs
            %   which minimizes overlap of t with other elements in ax. target is a two
            %   element array containing normalized coordinates [tx, ty]. If multiple
            %   locations are able to minimize overlap, PLACETEXT chooses the location
            %   which places the center of t closest to target.
            %
            %   [norm, xco, yco] = PLACETEXT(as, string, target, fs) also returns the X
            %   and Y coordinates before they are normalized to the axis limits.
            %
            % Adapted from:
            %
            % Peter Mao (2020). textbp: text with legend-style "best" placement
            % (https://www.mathworks.com/matlabcentral/fileexchange/11466-textbp-text-with
            % -legend-style-best-placement), MATLAB Central File Exchange.
            
            if isempty(s)
                % If supplied string is empty, placement does not matter.
                norm = [0 0];
                xco = 0;
                yco = 0;
                return
            end
            
            try
                temp = text(ax, 0, 0, s,...
                    'Units', 'normalized',...
                    'FontSize', fs); % Temporarily draw s on the axes.
                w = temp.Extent(3); % Width of s, normalized.
                h = temp.Extent(4); % Height of s, normalized.
                delete(temp); % Clear the drawing.

                oglimx = ax.XLim;
                oglimy = ax.YLim;
                W = diff(oglimx); % Width of the axes.
                H = diff(oglimy); % Height of the axes.
                twdt = W*w; % Width of s, in data points.
                thgt = H*h; % Height of s, in data points.

                limx = oglimx + [0.03*W -0.03*W];
                limy = oglimy + [0.03*H -0.03*H];
                xcoord = limx(1):W/100:limx(2)-twdt; % Fairly large number of evenly spaced X coordinates...
                ycoord = limy(1):H/100:limy(2)-thgt; % ...and evenly spaced Y coordinates as candidates.

                if isempty(xcoord) || isempty(ycoord)
                    % String is too big to fit on axes, so exit and return -1.
                    norm = -1;
                    xco = -1;
                    yco = -1;
                    return
                end

                kids = ax.Children; % All elements that are currently on these axes.
                x = cell(1, numel(kids));
                y = cell(1, numel(kids));
                for i = 1:numel(kids)
                    % For each element...
                    kid = kids(i);

                    % ...determine the X and Y coordinates which (roughly) trace out the
                    % element. Only considers line plots and text elements.
                    xtemp = [];
                    ytemp = [];
                    if strcmp(kid.Type, 'line')
                        xtemp = kid.XData;
                        ytemp = kid.YData;
                    elseif strcmp(kid.Type, 'text')
                        tmpunits = kid.Units;
                        kid.Units = 'data';
                        ext = kid.Extent;
                        kid.Units = tmpunits;
                        N = ceil(ext(3)*100/W);
                        M = ceil(ext(4)*100/H);
                        xtemp = repmat(linspace(ext(1), ext(1)+ext(3), N), 1, M);
                        ytemp = reshape(repmat(linspace(ext(2), ext(2)+ext(4), M), N, 1), 1, []);
                    end
                    x{i} = xtemp(:)'; % Force x{i} to be row vector.
                    y{i} = ytemp(:)'; % Force y{i} to be row vector.
                end
                x = [x{~cellfun(@isempty, x)}];
                y = [y{~cellfun(@isempty, y)}];

                pop = zeros(numel(ycoord), numel(xcoord));
                j = 1;
                for yp = ycoord
                    i = 1;
                    for xp = xcoord
                        % For each candidate point (xp, yp), determine how many of the above
                        % points fall within the box that has left lower corner (xp, yp),
                        % width twdt, and height thgt.
                        pop(j, i) = sum(sum((x >= xp).*(x <= xp+twdt).*(y >= yp).*(y <= yp+thgt)));
                        i = i+1;
                    end
                    j = j+1;
                end

                % Logical array the same size as pop with true where pop is minimized and
                % false elsewhere.
                popmin = pop == min(min(pop));

                if sum(sum(popmin)) > 1
                    % If there are multiple minima within pop (i.e. multiple best choices):

                    % 'Blur' pop. For any point (i,j), add the values at (i-1,j), (i+1,j),
                    % (i,j-1), and (i,j+1).
                    [a, b] = size(pop);
                    if a == 1
                        % If there are no vertical neighbors...
                        pop = [pop(:,2),pop(:,1:(b-1))]+[pop(:,2:b),pop(:,(b-1))] + pop;
                    elseif b == 1
                        % If there are no horizontal neighbors...
                        pop = [pop(2,:)',pop(1:(a-1),:)']'+[pop(2:a,:)',pop((a-1),:)']' + pop;
                    else
                        % If there are neighbors in both dimensions...
                        pop = [pop(2,:)',pop(1:(a-1),:)']'+[pop(2:a,:)',pop((a-1),:)']'+...
                        [pop(:,2),pop(:,1:(b-1))]+[pop(:,2:b),pop(:,(b-1))] + pop;
                    end

                    % Logical array the same size as pop with true where the locations of the
                    % previous minima in pop are still minima in the newly blurred pop and
                    % false elsewhere.
                    popx = popmin.*(pop == min(pop(popmin)));

                    if sum(sum(popx)) > 1
                        % If there are multiple minima within the blurred pop:

                        % Find normalized X and Y coordinates of all minima.
                        [i, j] = find(popx);
                        x = (xcoord(j)-oglimx(1))/W;
                        y = (ycoord(i)-oglimy(1))/H;

                        try
                            % If target is valid, use it.
                            tx = target(1);
                            ty = target(2);
                        catch
                            % Otherwise, set target to (0, 0).
                            tx = 0;
                            ty = 0;
                        end

                        % Find the left lower coordinate of the box containing target.
                        tx = tx - twdt/(2*W);
                        ty = ty - thgt/(2*H);

                        % Find the minima of the blurred pop which is closest to target.
                        [~, m] = min((x-tx).^2 + (y-ty).^2);
                        i = i(m);
                        j = j(m);
                    else
                        % If there is only one minimum within the blurred pop, consider it the
                        % optimal location.
                        [i, j] = find(popx);
                    end
                else
                    % If there is only one minimum within pop, consider it the optimal
                    % location.
                    [i, j] = find(popmin);
                end

                % Convert indices to plot coordinates.
                xco = xcoord(j);
                yco = ycoord(i);

                % Normalize the coordinates against the axes limits.
                norm = [(xco-oglimx(1))/W, (yco-oglimy(1))/H];
            catch
                norm = -1;
                xco = -1;
                yco = -1;
            end
        end
        
        function varargout = dealCells
            % DEALCELLS Sets each output variable to a string containing that variable's
            % order among all of the output variables.
            
            varargout = compose('%d', 1:nargout)';
        end
        
        function [x, y, z] = findClosestPoint(p1, p2, X, Y, Z)
            % FINDCLOSESTPOINT Finds the point on a surface which is closest to the line
            % between two other points.
            %   [x, y, z] = FINDCLOSESTPOINT(p1, p2, X, Y, Z) returns the
            %   three-dimensional coordinates of a point P on a surface S which is defined
            %   by XData X, YData Y, and ZData Z such that P is closer to the line between
            %   p1 and p2 than any other point on S. If S is plotted on axes ax and p1 and
            %   p2 are the points where the mouse position enters and exits ax,
            %   respectively, then P gives the point on S which is closest to the user's
            %   mouse position.
            
            if ~isvector(X)
                X = X(1,:);
            end
            if ~isvector(Y)
                Y = Y(:,1);
            end
            d = p2 - p1;
            Z2 = ((X(:)).' - p1(1))/d(1) - (Y(:) - p1(2))/d(2);
            [~, ind] = min(abs(Z2), [], 2);
            N = size(Z2, 1);
            Z3 = Z(sub2ind(size(Z),ind',1:N));
            i = find(diff(sign((Z3-p1(3))/d(3) - (X(ind)-p1(1))/d(1))), 1);
            y = Y(i);
            x = X(ind(i));
            z = Z3(i);
        end
        
        function changeData(p, x, y, z)
            % CHANGEDATA Updates XData, YData, ZData of a plotted element.
            %   CHANGEDATA(p, p2) sets the XData, YData, and ZData of p to the XData,
            %   YData, and ZData of p2.
            %
            %   CHANGEDATA(p, x, y, z) sets the XData of p to x, the YData of p to y, and
            %   the ZData of p to z.
            
            if nargin == 2
                p.XData = x.XData;
                p.YData = x.YData;
                p.ZData = x.ZData;
                x.Visible = 'off';
            else
                p.XData = x;
                p.YData = y;
                if nargin == 4
                    p.ZData = z;
                end
            end
            p.Visible = 'on';
        end
        
        function [in, p1, p2] = inAxes(ax, ex)
            % INAXES Determines whether the mouse is within the limits of a set of axes.
            %   in = INAXES(ax) checks each axes in ax and returns an array with elements
            %   set to true where the mouse is within the limits of the axes and elements
            %   set to false where the mouse is not within the limits of the axes.
            %
            %   in = INAXES(ax, ex) allows for the X coordinate of the mouse to be ex data
            %   points outside of the X limits of each axes in ax.
            %
            %   [in, p1, p2] = INAXES(ax, ___ ) also returns arrays p1 and p2 which, for
            %   each axes in ax, contain the coordinates where the mouse enters and exits
            %   the axes, respectively.
            
            e = 1e-6;
            if nargin < 2
                ex = e;
            end
            ey = e;
            ez = e;
            in = false(length(ax), 1);
            p1 = zeros(length(ax), 3);
            p2 = zeros(length(ax), 3);
            pt = cell(length(ax), 1);
            for i = 1:length(ax)
                pt{i} = ax(i).CurrentPoint;
            end
            for i = 1:length(ax)
                lim = axis(ax(i));
                if length(lim) ~= 6
                    lim(5:6) = ax(i).ZLim;
                end
                p = pt{i};
                p1(i,:) = p(1,:);
                p2(i,:) = p(2,:);
                in(i) = p1(i,1) >= lim(1)-ex && p1(i,1) <= lim(2)+ex &&...
                        p1(i,2) >= lim(3)-ey && p1(i,2) <= lim(4)+ey &&...
                        p1(i,3) >= lim(5)-ez && p1(i,3) <= lim(6)+ez &&...
                        p2(i,1) >= lim(1)-ex && p2(i,1) <= lim(2)+ex &&...
                        p2(i,2) >= lim(3)-ey && p2(i,2) <= lim(4)+ey &&...
                        p2(i,3) >= lim(5)-ez && p2(i,3) <= lim(6)+ez;
            end
        end
        
        function [col, i] = getExcelCol(i)
            % GETEXCELCOL Returns the label of an excel column.
            %   col = GETEXCELCOL(i) returns a string col with the same length as i where
            %   each letter in col is the ith letter of the alphabet.
            %   col = GETEXCELCOL([1 26 26]), for example, returns col = 'AZZ'.
            %
            %   [col, i] = GETEXCELCOL(i) also returns i, the input to GETEXCELCOL which
            %   would give the label of the next adjacent Excel column. Excel columns are
            %   labeled alphabetically where the column after 'Z' is 'AA' and the column
            %   after 'AZZ' is 'BAA'. [col, i] = GETEXCELCOL([1 26 26]), then, returns
            %   col = 'AZZ' and i = [2 1 1].
            
            letters = 'A':'Z';
            col = letters(i);
            if nargout == 2
                i(end) = i(end) + 1;
                while any(i > length(letters))
                    j = find(i > length(letters), 1, 'last');
                    if j == 1
                        i = [ones(size(i)) 1];
                    else
                        i(j) = 1;
                        i(j-1) = i(j-1)+1;
                    end
                end
            end
        end
        
        function output = closestTrue(input)
            % CLOSESTTRUE Returns the index of the closest true value to each position.
            %   output = CLOSESTTRUE(input) returns a vector output with the same length
            %   as input, N, such that for any i in [1, N], output(i) is the closest value
            %   to i for which input(output(i)) is true.
            
            input = input(:).';
            i = find(input);
            [~, j] = min(abs((1:length(input)) - i'), [], 1);
            output = i(j);
        end
        
        function resetAxesDefaultLimits(ax)
            % Removes the 'matlab_graphics_resetplotview' key and corresponding value from
            % the appdata of ax. This key's value seems to be set once the user first
            % interacts with the axes (e.g. through panning) and determines the axes
            % limits when 'Reset View' is clicked. If we programmatically change the axes
            % limits after the user has interacted with the axes, we want to reset these
            % default axes limits so that clicking 'Reset View' does not take the axes
            % back to their previous limits.
            
            try
                rmappdata(ax, 'matlab_graphics_resetplotview');
            catch
            end
        end
        
        function n = numSheets(file)
            % NUMSHEETS Returns number of sheets in an Excel workbook.
            %   n = NUMSHEETS(file) returns the number of sheets contained within the
            %   Excel workbook specified by the character array file.
            
            try
                % 2019b and later:
                sheets = sheetnames(file);
                n = length(sheets);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                    % 2019a and earlier:
                    if ispc
                        % Windows
                        [~, sheets] = xlsfinfo(file);
                        n = length(sheets);
                    else
                        % Mac (and probably also Linux?)
                        warning('off', 'MATLAB:xlswrite:AddSheet')
                        lastwarn('');
                        s = 1;
                        [~, name] = fileparts(file);
                        if ~strcmpi(name, 'temp')
                            copyName = 'temp.xlsx';
                        else
                            copyName = 'temp_.xlsx';
                        end
                        copyfile(file, copyName);
                        while 1
                            writetable(table(""), copyName, 'Sheet', s, 'WriteVariableNames', 0);
                            [~, id] = lastwarn;
                            if strcmp(id, 'MATLAB:xlswrite:AddSheet')
                                delete(copyName);
                                n = s-1;
                                warning('on', 'MATLAB:xlswrite:AddSheet')
                                return
                            end
                            s = s+1;
                        end
                    end
                else
                    rethrow(ME)
                end
            end
        end
        
        function W = fillTable(T, D, mH)
            % FILLTABLE Fills a table while optimizing column widths.
            %   W = FILLTABLE(T) adjusts the columns of table T such that the sum of
            %   column widths equals the width of the table and, if possible, such that no
            %   column is too thin for the data it contains. Any extra width is
            %   distributed proportionally among all columns. The final column widths are
            %   returned in W.
            %
            %   W = FILLTABLE(T, D) uses data D to fill the table T, rather than using the
            %   data already present in the table.
            %
            %   W = FILLTABLE(T, D, mH) accepts mH, the maximum height of the table. If
            %   the data within D contains enough rows to exceed mH, a scrollbar will be
            %   added. This scrollbar affects the column widths and is assumed to have a
            %   width of 16 pixels.
            
            % Set table data.
            if nargin < 2
                % If user is not supplying new data, use data that is already present.
                D = T.Data;
            end
            if nargin < 3
                mH = Inf;
            end
            T.Data = D;
            
            % Initialize temporary text and table elements for measuring various sizes.
            tempText = uicontrol('Style', 'text',...
                'FontSize', T.FontSize,...
                'FontName', T.FontName,...
                'Visible', 'off');
            tempTable = uitable('Visible', 'off');
            
            % Determine minimum widths required to display all data in each column.
            minWidths = zeros(size(D,2),1);
            for c = 1:length(minWidths)
                % For each column...
                for r = 1:length(D(:,c))
                    % ...for each cell in that column...
                    if ischar(D{r,c})
                        % ...if the data in that cell is a character array...
                        tempText.String = D{r,c};
                        w = tempText.Extent(3); % ...determine it's width...
                        if w > minWidths(c)
                            % ...and update minWidths if necessary...
                            minWidths(c) = w;
                        end
                    end
                end
                % ...also consider the table column names:
                tempTable.ColumnName = T.ColumnName{c};
                if tempTable.Extent(3) > minWidths(c)
                    minWidths(c) = tempTable.Extent(3);
                end
            end
            
            % Update column height and widths.
            hgt = min([T.Extent(4) mH]);
            if hgt == mH
                % Table was too big, so slider was added.
                addSlider = true;
            else
                addSlider = false;
            end
            pUnits = T.Parent.Units;
            T.Parent.Units = 'pixels';
            hgt = hgt + 1/T.Parent.Position(4);
            T.Parent.Units = pUnits;
            T.Position([2 4]) = [T.Position(2)+T.Position(4)-hgt hgt];
            origUnits = T.Units;
            T.Units = 'pixels';
            if addSlider
                wdt = T.Position(3) - 16;
            else
                wdt = T.Position(3);
            end
            if sum(minWidths)+size(D,2)+1 <= wdt
                % If the sum of the minimum widths do not exceed the table width,
                % distribute the extra width proportionally among all columns.
                minWidths = minWidths + (wdt - (sum(minWidths)+size(D,2)+1)).*(minWidths./sum(minWidths));
                if sum(floor(minWidths))+size(D,2)+1 <wdt
                    % Final adjustments due to rounding.
                    N = floor(wdt - (sum(floor(minWidths))+size(D,2)+1));
                    minWidths(1:N) = minWidths(1:N)+1;
                end
                widths = num2cell(minWidths)';
                
                % Center all cells which contain character arrays, using html.
                Dc = cell(size(D));
                for c = 1:size(D,2)
                    if ischar(D{1,c})
                        Dc(:,c) = strcat(sprintf('<html><tr align=center><td width=%d>', widths{c}), D(:,c));
                    else
                        Dc(:,c) = D(:,c);
                    end
                end
                
                % Update the table data and column widths.
                T.Data = Dc;
                T.ColumnWidth = widths;
            else
                % Otherwise, the sum of the minimum widths required to display all data is
                % too large for the table, so resort to automatic width calculation.
                T.ColumnWidth = 'auto';
            end
            W = T.ColumnWidth;
            
            % Clean up.
            delete(tempText);
            delete(tempTable);
            T.Units = origUnits;
            
            T.Visible = 'on';
        end
    end
end
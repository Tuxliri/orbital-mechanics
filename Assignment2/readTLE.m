function [DATES, KEP]=readTLE(fileNAME, print)
% readTLE computes satellite ephemeris data from a list of two-line element (TLE) file, extracted from space-track.org database.
%
% PROTOTYPE:  [DATES, KEP]=readTLE(fileNAME, print)
%
% INPUTS:
%   • fileNAME: path to .TLE file (any standard two-line element file)
%                       NOTE: when download the TLE list copy the list in a.TXT file and change the extension to .TLE.
%   • (OPTIONAL) print: 1 or 0 to print on the command window the values extracted
%
% OUTPUTS: 
%   • DATES [Nx6]: date vector [yyyt, mm, dd, hh, mm, ss], N depends on the length of the list
%   • KEP [Nx7]: keplerian elements
%         KEP(k,1)=a : semimajor axis , km
%         KEP(k,2)=e : eccentricity , - 
%         KEP(k, 3)=i : inclination, rad
%         KEP(k, 4)=OMEGA : RAAN, rad
%         KEP(k, 5)=omega: argument of the perigee, rad
%         KEP(k,6)=theta: true anomaly, rad
%         KEP(k,7)=M;
% CALLED FUNCTIONS:
%   • result = chksum(str) : declared in this file
%
% CONTRIBUTORS
%   • Andrea Sportillo (andrea.sportillo@mail.polimi.it - Politecnico di Milano - Orbital Mechanics Course AA 2020/2021)
%      Based on script of Brett Pantalone from  North Carolina State
%      University: https://fr.mathworks.com/matlabcentral/fileexchange/56904-read-two-line-element-ephemeris-files
%
% LOG:
%   • 02/01/2021 v0: first issue
%       TBD : check calculation and initialize output variables to spped up the code
% NOTE: PLEASE, if you identify any bug in the variable calculation report it AS SOON AS POSSIBLE by email.

% OPTIONAL : print 
  if nargin<2
      print=0; % No print y default
  end
  
% OPEN FILE
  fd = fopen(fileNAME,'r');
  
  if fd < 0
      fd = fopen([fileNAME '.tle'],'r');
  end
  
  assert(fd > 0,['Can''t open file ' fileNAME ' for reading.']); % IN case of problems with the file

  A1 = fgetl(fd); % TLE 1
  A2 = fgetl(fd); % TLE 2
  
  k=0; % init
  
   SATNUM1 = str2double(A1(3:7)); % object id number
   
   % MAIN CYCLE TO EXTRACT DATA
  while ischar(A2) % until the end of the list in the file
    k=k+1;
      satnum = str2double(A1(3:7)); % object id number
      
    if satnum~= SATNUM1
        warning('The Object numer has changed, in this file there are TLE for more than 1 object!')
    end
    
    % Recover the date
        epoch(k)=str2double(A1(19:32)); % [YYDDD.fff]
        yy= str2double(A1(19:20)); % year
        yDay=str2double(A1(21:32));
         if yy>56  % no mission before sputnik 1 (1957)
            year=1900+yy;
            date=datevec(datetime([year 01 01 0 0 0])+caldays(fix(yDay)-1)+days(yDay-fix(yDay)));
        else
            year=2000+yy;
            date=datevec(datetime([year 01 01 0 0 0])+caldays(fix(yDay)-1)+days(yDay-fix(yDay)));
         end
        % Recover the keplerian elements
            n = str2double(A2(53:63));
            T = 86400/n; % Period, s
            a =1e-3* ((T/(2*pi))^2*398.6e12)^(1/3);
            e = str2double(['.' A2(27:33)]); % in TLEs orbits are assumed to be closed (e<1)
            i = str2double(A2(9:16))*pi/180;
            OMEGA = str2double(A2(18:25))*pi/180;
            omega = str2double(A2(35:42))*pi/180;
            M = str2double(A2(44:51))*pi/180;      
      
           % Recover the True anomaly from the mean Anomaly trhough the
           % Kepler's equation (Newton's Method)
            E0=M+( e*sin(M) ) / (1-sin(M+e)+sin(M) );  % rad 
            itMax=50;  tol=1e-6;
            it=0; err=1;
            E=E0;
        while err>tol && it<itMax
            it=it+1;
            Eold=E;
            E=Eold-(M-Eold+e*sin(Eold))/(e*cos(Eold)-1);  % x1=x0+f(x0)/f'(x0)
            err= abs(E-Eold)/abs(Eold);
        end
        
         % Compute the correct angular osition in the plane:
                 cosTheta=(cos(E)-e)/(1-e*cos(E)); 
                 sinTheta=sqrt(1-cosTheta^2); % Just the positive valueof sin(Theta), we'll check the sign with sin(M(i)) later on

                if  sin(M)>=0
                        theta=atan2(sinTheta, cosTheta);
                else
                        theta=atan2(-sinTheta, cosTheta);   
                end
                
         % OPTIONAL: print the current extracted values on the Cmatlabs
         % Command Window
      if print
          fprintf('%s\n', repmat('-',1,50));
          fprintf('Catalog Number: %d\n', satnum)
          fprintf('Epoch time: %s\n', A1(19:32)) 
          fprintf('Date: %s\n', num2str(date, 4) )
          fprintf('Semi-major axis: %.0f km\n', a)
          fprintf('Eccentricity: %f\n', e)
          fprintf('Inclination: %f rad\n', i)
          fprintf('RA of ascending node: %f rad\n', OMEGA)
          fprintf('Arg of perigee: %f rad\n', omega)
          fprintf('True anomaly: %f rad\n', theta)
          fprintf('Mean anomaly: %f rad\n', M)
          fprintf('Mean motion: %f rev/day\n', n)
          fprintf('Period of rev: %.0f s/rev\n', T)
      end
      
      % Check the 2 lines:
           assert(chksum(A1), 'Checksum failure on line 1');  % if error in the TLE
           assert(chksum(A2), 'Checksum failure on line 2');
           
        % Save Data for output
        DATES(k, :)=date;
        KEP(k,1)=a;
        KEP(k,2)=e;
        KEP(k, 3)=i;
        KEP(k, 4)=OMEGA;
        KEP(k, 5)=omega;
        KEP(k,6)=theta;
        KEP(k,7)=M;
        
        % Next 2 lines
        A1 = fgetl(fd);
        A2 = fgetl(fd);
        
  end
  
  % Erase duplicated lines
    [~,uidx] = unique(epoch,'stable'); % To erase duplicated rows (same day value)
    DATES= DATES(uidx,:);
    KEP=KEP(uidx,:);
  
  fclose(fd); % close file
end

%%
% Checksum (Modulo 10)
% Letters, blanks, periods, plus signs = 0; minus signs = 1
function result = chksum(str)
  result = false; c = 0;
  
  for k = 1:68
    if str(k) > '0' && str(k) <= '9'
      c = c + str(k) - 48;
    elseif str(k) == '-'
      c = c + 1;
    end
  end

  if mod(c,10) == str(69) - 48
    result = true;
  end
  
end

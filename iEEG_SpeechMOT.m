function targdist_fmri_final
%cd('C:\Users\Axolotl\Desktop\iEEG')
%Screen('Preference','SkipSyncTests', 1); %keep during debugging
%% written by Jason Rajsic, March 2011; Harrison Ritz, March 2015 (harrison.ritz@gmail.com); Bruno Mesquita, December 2022
%% Program requires participants to track 1:4 objects while listening to speech



%% Set-up Digital Triggers

% Get Device Info (Won't be used after we have the proper DeviceID)
dev = daqlist("ni");
DevInf = dev.DeviceInfo(1);

% Creting DataAcquisition Object
d = daq("ni")

% Adding Output Channels (The device name will depend on how it is set up)

addoutput(d,"Dev1","mot","Digital");
addoutput(d,"Dev1","sent","Digital");



%% set variables for key codes
KbName('UnifyKeyNames');
kbtrig = KbName('t');
kbstart = KbName('0');

kb1 = KbName('b');
kb2 = KbName('y');
kb3 = KbName('g');



%% get participant info
valid_file = 0;
while valid_file == 0,
    
    prompt = {'Participant #', 'CB (1:6)', 'Age', 'Sex'};
    answer = inputdlg(prompt, 'Caesar Set-up',1,{'1','1','23','F'});
    participantNumber = str2double(answer{1});
    cb_ans = str2double(answer{2});
    age = str2double(answer{3});
    sex = answer{4};
    
    file = sprintf('results_testing/targdist_MOT_%d.csv', participantNumber);
    
    
    valid_file = 1;
    if fopen(file) ~= -1,
        choice = questdlg('You are about to overwrite an existing data file. Are you sure you want to do this?','Error','Yes','No','No');
        if strcmp(choice, 'No') == 1,
            valid_file = 0;
        end
    end
end


[resultsFile, ~] = fopen(file,'r+');
% fprintf(resultsFile, 'Caesar Results');
% fprintf(resultsFile, '\n\n,,Attend:,0 = speech,1 = MOT');
% fprintf(resultsFile, '\n,,Speech:,0 = clear, 1 = NV6, 2 = NV12');

expOnset = GetSecs;
tic;

fprintf(resultsFile,'pt,Age,Sex,CB,Block,Trial,cond,attend,speech,sent num,Sentence,mot file,targets,num Dots,exp onset,block trigger onset,block onset,trial que onset,movie onset,query onset,ITI onset,trigger,trigger time,Correct Response,key press,key time,acc,rt\n'); % make a header row for subject file
%PsychDebugWindowConfiguration(0,0.5);
Screen('Preference','SkipSyncTests', 1);
 % Find out how many screens and use largest screen number.
    whichScreen = 1;

    % Open a new window.
[window,rect] = Screen('OpenWindow', whichScreen);
rng('shuffle'); % resets the random number generator
pause on;
%HideCursor; %disable during debugging

%% Experiment Parameters
nBlocks = 6;        %number of blocks
nTrials = 24;       %number of trials/block
cb_len = 4;         %number of task conditons
cb_tr = 36;         %number of trials for each task x speech (block*trial / cb_len)
breakBlock = 3;     %number of blocks until break
blankNum = 0;      %blank trials
rotNum = 16;        %dummy trials

probeNum = 3; %number of probes

%timing parameters
fixlen =    0.3;      % ITI
targetDur = 1.8;      % duration of target cues
trackDur =  5.0;      % duration of tracking in seconds
pauseDur =  0.1;      % time between end of movie and probe
respDur =   2.8;      % time to respond to probe
trigOffset = respDur - (4 - (trackDur/2 + pauseDur)); %duration of the response window after an optimal (4 sec after sentence midpoint) trigger 



% display parameters
trackSize = [700,700];  % size of tracking space
objSize = 14;           %size, in pixels, of dot
dist = 12;              % number distractors
frameRate = 60;     %frames per second
fixSize = 4;        %size of the fixation box (well, half of it - the distance from the middle to edge horizontally/vertically)
qFrame = 1;         %width of frame around each tracking quadrant
%Fs = 44100;
%shrink = 1;         % reduce the real size of the movie by some amount to conserve memory.
%screen stuff
rect = Screen('Rect',window)/1;
xMid = rect(3)/2;
yMid = rect(4)/2;
white = [255,255,255];
black=[0,0,0];
grey = [128,128,128];
green = [0,255,0];
red = [255,0,0];
blue = [0,128,255];
yellow = [255,255,0];
backgroundColour = [100,100,100];
Screen('TextFont', window, 'Arial');
fixCol = yellow;
frameCol = grey;
trackRect = [xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2];


kb = vertcat([kb1,1],[kb2,2],[kb3,3]);
accuracy = nan(nBlocks,nTrials);

%% set condition: cb(:,2) = attend (0 = speech, 1 = mot); cb(:,3) = speech (0 = clear, 1 = nv12)
cb = [[0,1,0,1]',[0,0,1,1]'];
switch cb_ans
    case 1
        cb = [[1,2,3,4]',cb];
    case 2
        cb = [[2,3,4,1]',cb];
    case 3
        cb = [[3,4,1,2]',cb];
    case 4
        cb = [[4,1,2,3]',cb];
end

%get datasource
[cond_list] = xlsread('conditions_file/targ_datasource_bm','cond','A2:C25'); %get CB info
[num, txt] = xlsread('conditions_file/targ_datasource_bm','set','A2:G145'); %get stim info

%randomize the order of the MOT movies
stim = randperm(nTrials*nBlocks);
count = 1;

%randomize order of sentences within each condtion
rand_sent = [];
for r = 1:cb_len
    rand_sent = [rand_sent, randsample((r-1)*cb_tr+1:r*cb_tr,cb_tr)'];
end
sentcount = ones(1,cb_len);

%% Set up BASELINE trials (dummy == baseline)
dummyNum = [2,3,3]; %number of each baseline type per block
%blankBl = [dummyNum(randperm(3)),dummyNum(randperm(3)),dummyNum(randperm(3))];
rotBl =   [dummyNum(randperm(3)),dummyNum(randperm(3))];

dotNum = [1,4]; %number of dots in the rotated baseline trials
dotOrd = dotNum(randperm(2));
for ii = 1:7
    dotOrd = [dotOrd,dotNum(randperm(2))];
end
dotCount = 1;

dSent = randperm(rotNum);
dSentCount = 1;


%% flip correction (singal flip early in the frame refresh, so will flip as soon
%as the monitor refreshes)
er =  1/61;

% set up initial audioplayer
load handel;
handel_test = audioplayer(y, Fs);
play(handel_test,[1, 10])
clear handel_test
Fs = 44100;


%% Start Screen
Screen('FillRect',window,backgroundColour);

inst1 = 'You will see several dots and be asked to either TRACK or LISTEN';
inst1W=TextBounds(window,inst1);
inst2 = 'If TRACK, focus on the red dots and select the target';
inst2W=TextBounds(window,inst2);
inst3 = 'If LISTEN, focus on the sentence and respond if you understood it';
inst3W=TextBounds(window,inst3);
inst4 = 'Press any key to start';
inst4W=TextBounds(window,inst4);

Screen('DrawText',window,inst1,xMid-(inst1W(3)-inst1W(1))/2,yMid-100,white);
Screen('DrawText',window,inst2,xMid-(inst2W(3)-inst2W(1))/2,yMid-50,white);
Screen('DrawText',window,inst3,xMid-(inst3W(3)-inst3W(1))/2,yMid-0,white);
Screen('DrawText',window,inst4,xMid-(inst4W(3)-inst4W(1))/2,yMid+50,white);

Screen('Flip',window);
KbWait([], 3);

%% BLOCKS
disp(nBlocks)
for block=1:nBlocks
    
       
    %blank screen while waiting for
    Screen('FillRect',window,backgroundColour);
    Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
    Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
    
    Screen('Flip',window); %present fixation
    
    %get first trigger
    keyCode = zeros(1,256);
    
    buttonTime = GetSecs; % check for time at start of expt
    
    
    blTrigOnset = buttonTime; %first trigger of the block
    blOnset = blTrigOnset + trigOffset + fixlen; %first trial onset
    trigSet = blOnset;
    adjBlOnset = blOnset; %reset adjusted block onset
    
    % Indicate block will begin
    Screen('FillRect',window,backgroundColour);
    Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
    Screen('FillRect',window,green,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
    Screen('Flip',window); %present start-of-block fixation
    
    %set response and order
    resp = randsample(vertcat(ones(nTrials/probeNum,1),ones(nTrials/probeNum,1)*2,ones(nTrials/probeNum,1)*3),nTrials); %randomize response order (valid resp: 1:3)
    order = horzcat(cond_list(randperm(nTrials),:),resp);
    
   
    %number of dummy trials this block
    %blankTrNum = blankBl(block);
    rotTrNum = rotBl(block);
    
    %determine which trials will procede blank and rot trials
    dummyTr = randperm(nTrials, rotTrNum);
    %blankTr = dummyTr(1 : blankTrNum); %blank trials
    rotTr = dummyTr(end-rotTrNum+1 : end); %rotated speech trials
    
    
    %% TRIALS
    for trial=1:nTrials
        
        
        %trial parameters
        cond = order(trial,2);
        attend = cb(cb(:,1)==cond,2);       %set attend (0 = speech, 1 = mot)
        speech = cb(cb(:,1)==cond,3);       %set speech (0= clear, 1 = nv12)
        numTrack = order(trial,3);          %set number of targets        
        corr_resp = order(trial,4);

        numDots = numTrack+dist;        %number of dots in total
        
        targets = randsample(numDots,numTrack);                         %pick targets from the dots
        corr = datasample(targets,1);                                   %pick a target to be queried
        foils = randsample(setxor(1:numDots,targets), probeNum-1);      %pick dots that are not targets to be foils
        
        
        %load clear or nv speech (I saved the y of a wav into a mat
        %file for faster loading.. its at the trial level so it doesnt
        %matter too much)        
        
        switch speech
            case 0 %file path
                load(fullfile('stimuli/cl_stim',sprintf('cl_%s',txt{rand_sent(sentcount(cond),cond),2})))
                
            case 1
                load(fullfile('stimuli/nv12_stim',sprintf('nv12_%s',txt{rand_sent(sentcount(cond),cond),2})))     
        end
        
        
        sentcount(cond) = sentcount(cond)+1;
        
        sent = audioplayer(y,Fs);
        clear y
        
        % load mot object locations
        % I make the movies with MOTmovie (via
        % MOTmaker), and it will give me a 4D matrix: frame by object by XY
        % coordinates. I save this list of locs to a mat file for faster loading.
        % Again, its at a trial level (otherwise you'd make the
        % location matrix at the beginning of each trial), but its cleaner and it easily lets
        % you see the locations of the objects for item effects (heaven
        % forbid). I made the movies 60fps, but I can use any
        % fps below that, it just wont make it through the entire movie
        
        
        %% mot file
        
        % motfile = sprintf('sp1_mot_targ_%.3d',stim(count));
        %motfile = sprintf('dt18_sp1_sz24_mot_%d',stim(count));
        %motfile = sprintf('dt18_bf32_mot_%d',stim(count));
        %motfile = sprintf('dt16_mot_%d',stim(count));
        %motfile = sprintf('sm_dt18_mot_%d',stim(count));
        motfile = sprintf('mot_%.3d',stim(count));
        load(fullfile('stimuli/mots_targ', motfile));
        count = count+1;
        
        %% MOVIE
        %% Present Objects  
        
        Screen('FillRect',window,backgroundColour); %% drawing the target cues
        Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
        
        temp = ones(4, numDots); %pre-allocate temp
        for o=1:numDots% make a matrix of object locations, columns are different objects, rows are each object's xy coords
            temp(:,o)= [locs(1,1,o,1)+trackRect(1)-objSize; locs(1,2,o,1)+trackRect(2)-objSize; locs(1,1,o,1)+trackRect(1)+objSize; locs(1,2,o,1)+trackRect(2)+objSize];
        end
        Screen('FillOval',window,white,temp,objSize); %draw the objects
        
        Screen('FrameRect',window,frameCol,trackRect,qFrame);
        
        %% draw targets
        
        if attend
            temp = ones(4,length(targets));
            for t=1:length(targets)
                temp(:,t) = [locs(1,1,targets(t),1)+trackRect(1)-objSize; locs(1,2,targets(t),1)+trackRect(2)-objSize; locs(1,1,targets(t),1)+trackRect(1)+objSize; locs(1,2,targets(t),1)+trackRect(2)+objSize];
            end
            Screen('FillOval',window,red,temp,objSize);
        end
        
        % It seems to be preferable to make a matrix of object
        % locations and drawing them all at once rather than looping
        % through drawing each oval
        
        
        Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
        
        %draw TRACK or LISTEN
        Screen('TextSize', window, 32);
        if attend
            Screen('DrawText', window, 'TRACK',xMid-73, yMid-52, blue);
        else
            Screen('DrawText', window, 'LISTEN',xMid-78, yMid+6, blue);
        end
        
        % trOnset = Screen('Flip',window,fixOn+fixlen-er);
        trOnset = Screen('Flip',window, adjBlOnset + ((trial-1)*10) - er);        
        %trOnset = Screen('Flip',window, trigSet - er);
        
        Screen('FillRect',window,backgroundColour);
        Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
        
        temp = [locs(1,1,o,1)+trackRect(1)-objSize; locs(1,2,o,1)+trackRect(2)-objSize; locs(1,1,o,1)+trackRect(1)+objSize; locs(1,2,o,1)+trackRect(2)+objSize];
        Screen('FillOval',window,white,temp,objSize);
        
        Screen('FrameRect',window,frameCol,trackRect,qFrame);
        
        % Play Sentence
        WaitSecs('UntilTime', trOnset + targetDur - (er*2));
        play(sent);
        
        % SEND SENT TRIGGER
        
        write(d,[2 0]); % IS THAT IT? NOT SURE
        
        movOnset = Screen('Flip',window, trOnset + targetDur - er);
        
        %ti = GetSecs;
        
        %% tracking movie begin
        
        % SEND MOT TRIGGER
        
        write(d,[1 0]); % IS THAT IT? NOT SURE
        
        %mov=tic;
        for frame = 1:frameRate*trackDur
            
            Screen('FillRect',window,backgroundColour);
            Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
            
            temp = ones(4, numDots);
            for o=1:numDots
                temp(:,o)= [locs(1,1,o,frame)+trackRect(1)-objSize; locs(1,2,o,frame)+trackRect(2)-objSize; locs(1,1,o,frame)+trackRect(1)+objSize; locs(1,2,o,frame)+trackRect(2)+objSize];
            end
            Screen('FillOval',window,white,temp,objSize);
            
            Screen('FrameRect',window,frameCol,trackRect,qFrame);
            Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
            
            
            %WaitSecs('UntilTime', ti + 1/frameRate-er);
            Screen('Flip',window);
            %ti = GetSecs;
            
            
        end
        
        
        %% display probe
        Screen('FillRect',window,backgroundColour);
        Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
        
        temp = ones(4, numDots);
        for o=1:numDots
            temp(:,o)= [locs(1,1,o,frame)+trackRect(1)-objSize; locs(1,2,o,frame)+trackRect(2)-objSize; locs(1,1,o,frame)+trackRect(1)+objSize; locs(1,2,o,frame)+trackRect(2)+objSize];
        end
        Screen('FillOval',window,white,temp,objSize);
        
        
        
        Screen('FrameRect',window,frameCol,trackRect,qFrame);
        Screen('TextSize', window, 12);
        
        % if attend mot, query probes (one is target), otherwise query gist
        if attend
            f = 1; %keep track of how many foils
            for p = 1:probeNum
                if p == corr_resp
                    
                    temp = [locs(1,1,corr,frame)+trackRect(1)-objSize; locs(1,2,corr,frame)+trackRect(2)-objSize; locs(1,1,corr,frame)+trackRect(1)+objSize; locs(1,2,corr,frame)+trackRect(2)+objSize];
                    Screen('FillOval',window,blue,temp,objSize);
                    Screen('DrawText',window,num2str(p),temp(1,1)+objSize-4, temp(2,1)+objSize-9,white);
                    
                    
                else
                    
                    temp = [locs(1,1,foils(f),frame)+trackRect(1)-objSize; locs(1,2,foils(f),frame)+trackRect(2)-objSize; locs(1,1,foils(f),frame)+trackRect(1)+objSize; locs(1,2,foils(f),frame)+trackRect(2)+objSize];
                    Screen('FillOval',window,blue,temp,objSize);
                    Screen('DrawText',window,num2str(p),temp(1,1)+objSize-4,temp(2,1)+objSize-9, white);
                    f = f+1;
                    
                    
                end
            end
            Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
            
        else
            Screen('TextSize', window, 32);
            
            Screen('DrawText', window,'GIST?',xMid-60, yMid-23, blue);
            Screen('DrawText', window,'YES',xMid/2 -42, yMid+150, yellow);
            Screen('DrawText', window,'NO',xMid*1.5 - 33, yMid+150, yellow);
        end
        
        %WaitSecs('UntilTime', ti + pauseDur-er);
        %startTime = Screen('Flip',window,ti + pauseDur - er);
        
        queryOnset = Screen('Flip',window, trOnset + targetDur + trackDur + pauseDur - er);
        
        %% get response
        
        % load ITI screen
        Screen('FillRect',window,backgroundColour);
        Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
        Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
        
        %reset response variables
        respcount = 1;
        key = zeros(1,256);
        time = [0];
        
        while buttonTime <= queryOnset + respDur - (er*2)
            
            [press,buttonTime,keyCode] = KbCheck(-1); % check for a keypress
            
            if press %if button press, record button and time
                key(respcount,:) = keyCode;
                time(respcount) = buttonTime;
                respcount = respcount+1;
            end
            
        end
        
        %% Show ITI screen while file things
        
        fixOn = Screen('Flip',window, queryOnset + respDur - er); %present fixation
        
        
        %get unique keypresses and their onset
        [key,ia] = unique(key,'rows','stable');
        time = time(ia);
        
        % was there a trigger, and at what time
        trigRow = find(key(:,kbtrig),1);
        trig = sum(trigRow > 0);
        
        if trig
            trigTime = time(trigRow);
            key = setxor(key,key(trigRow,:),'rows','stable'); %remove trigger from keypresses
            time = setxor(time,time(trigRow));
            %trigSet = trigTime + trigOffset + fixlen; %first trial onset     
            
    
        else              
            trigTime = NaN;
            %trigSet = adjBlOnset + ((trial-1)*10); %first trial onset
        end
        
        %was there a response, when, and calc rt
        
        if ~sum(key)
            respOnset = NaN;
            RT = 0;
            key = 1;
            keyPress = 0;
        else
            respOnset = time(end);
            RT = respOnset - queryOnset;
            keyPress = kb(kb==find(key(end,:),1),2);
        end
        
        % if attending to MOT, match response to correct response
        % if attending to speech, record whether understood gist
        
        if attend %if MOT
            if find(key(end,:),1) == kb(corr_resp)
                accuracy(block,trial) = 1;
            else
                accuracy(block,trial) = 0;
            end
        else %if speeech
            corr_resp = 1;
            if find(key(end,:),1) == kb(1,1)
                accuracy(block,trial) = 1;
            else
                accuracy(block,trial) = 0;
            end
            
        end
        
        
        
        %print results      exp onset,block trigger onset,block onset,  trial/que onset,movie onset,query onset,  ITI onset,trigger,    trigger time,Correct Response,  key press,key time,acc,rt
        fprintf(resultsFile,'%d,%d,%s,%d,%d,%d,%d,%d,%d,%d,%s,%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%d,%d,%.4f,%d,%.4f\n', participantNumber,age,sex,cb_ans,block,trial,cond,attend,speech,num(rand_sent(sentcount(cond)-1,cond),1),txt{rand_sent(sentcount(cond)-1,cond),1},motfile,numTrack,numDots,expOnset,blTrigOnset,blOnset,trOnset,movOnset,queryOnset,fixOn,trig,trigTime,corr_resp,keyPress,respOnset,accuracy(block,trial),RT);
        clear locs sent
        
        
  %% Rotated Trials  ~~ Listen trials with rotated speech
        
        if max(rotTr == trial)
            
            
            
            %get number of dots
            numTrack = dotOrd(dotCount);          %set number of targets
           
            numDots = numTrack+dist;              %number of dots in total
            dotCount = dotCount + 1;
            
            %% load speech file
            load(fullfile('r_stim',sprintf('r_stim_%.2d',dSent(dSentCount))));
            dSentCount = dSentCount + 1;
            
            sent = audioplayer(y,Fs);
            clear y
            
            %% load mot file
            motfile = sprintf('d_mot_%.3d',numTrack);
            
            
            load(fullfile('stimuli/mots_targ', motfile));
            
            
            %% draw screen
            Screen('FillRect',window,backgroundColour); %% drawing the target cues
            Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
            
            temp = ones(4, numDots); %pre-allocate temp
            for o=1:numDots% make a matrix of object locations, columns are different objects, rows are each object's xy coords
                temp(:,o)= [locs(1,1,o,1)+trackRect(1)-objSize; locs(1,2,o,1)+trackRect(2)-objSize; locs(1,1,o,1)+trackRect(1)+objSize; locs(1,2,o,1)+trackRect(2)+objSize];
            end
            
            Screen('FillOval',window,white,temp,objSize); %draw the objects
            Screen('FrameRect',window,frameCol,trackRect,qFrame);
            
            %% draw LISTEN (?)
            Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
            Screen('TextSize', window, 32);
            Screen('DrawText', window, 'LISTEN',xMid-78, yMid+6, blue);
            
            trOnset = Screen('Flip',window, adjBlOnset + (trial*10) - er);
            %trOnset = Screen('Flip',window, trigSet - er);
            
            Screen('FillRect',window,backgroundColour);
            Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
            
            temp = [locs(1,1,o,1)+trackRect(1)-objSize; locs(1,2,o,1)+trackRect(2)-objSize; locs(1,1,o,1)+trackRect(1)+objSize; locs(1,2,o,1)+trackRect(2)+objSize];
            Screen('FillOval',window,white,temp,objSize);
            
            Screen('FrameRect',window,frameCol,trackRect,qFrame);
            
            %% play sentence
            WaitSecs('UntilTime', trOnset + targetDur - (er*2));
            play(sent);
            
            % SEND SENT TRIGGER
        
        write(d,[2 0]); % IS THAT IT? NOT SURE
        
            movOnset = Screen('Flip',window, trOnset + targetDur - er);
            
            %% movie
            
            % SEND MOT TRIGGER
        
        write(d,[1 0]); % IS THAT IT? NOT SURE
        
        
            for frame = 1:frameRate*trackDur
                
                Screen('FillRect',window,backgroundColour);
                Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
                
                temp = ones(4, numDots);
                for o=1:numDots
                    temp(:,o)= [locs(1,1,o,frame)+trackRect(1)-objSize; locs(1,2,o,frame)+trackRect(2)-objSize; locs(1,1,o,frame)+trackRect(1)+objSize; locs(1,2,o,frame)+trackRect(2)+objSize];
                end
                Screen('FillOval',window,white,temp,objSize);
                
                Screen('FrameRect',window,frameCol,trackRect,qFrame);
                Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
                
                Screen('Flip',window);
                
            end
            
            
            %% display probes
            Screen('FillRect',window,backgroundColour);
            Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
            
            temp = ones(4, numDots);
            for o=1:numDots
                temp(:,o)= [locs(1,1,o,frame)+trackRect(1)-objSize; locs(1,2,o,frame)+trackRect(2)-objSize; locs(1,1,o,frame)+trackRect(1)+objSize; locs(1,2,o,frame)+trackRect(2)+objSize];
            end
            Screen('FillOval',window,white,temp,objSize);
                        
            Screen('FrameRect',window,frameCol,trackRect,qFrame);
            Screen('TextSize', window, 12);
            
            %% query gist
            
            Screen('TextSize', window, 32);
            
            Screen('DrawText', window,'GIST?',xMid-60, yMid-23, blue);
            Screen('DrawText', window,'YES',xMid/2 -42, yMid+150, yellow);
            Screen('DrawText', window,'NO',xMid*1.5 - 33, yMid+150, yellow);
            
            
            queryOnset = Screen('Flip',window, trOnset + targetDur + trackDur + pauseDur - er);
            
            %% get response
            
            % load ITI screen
            Screen('FillRect',window,backgroundColour);
            Screen('FillRect',window,black,[xMid-trackSize(1)/2,yMid-trackSize(2)/2,xMid+trackSize(1)/2,yMid+trackSize(2)/2]);
            Screen('FillRect',window,fixCol,[xMid-fixSize,yMid-fixSize,xMid+fixSize,yMid+fixSize]); % draw fixation box
            
            %reset response variables
            respcount = 1;
            key = zeros(1,256);
            time = [0];
            
            while buttonTime <= queryOnset + respDur - (er*2)
                [press,buttonTime,keyCode] = KbCheck; % check for a keypress
                
                if press %if button press, record button and time
                    key(respcount,:) = keyCode;
                    time(respcount) = buttonTime;
                    respcount = respcount+1;
                end
                
            end
            
            %% Show ITI screen while file things
            fixOn = Screen('Flip',window, queryOnset + respDur - er); %present fixation
            
            %get unique keypresses and their onset
            [key,ia] = unique(key,'rows','stable');
            time = time(ia);
            
            % was there a trigger, and at what time
            trigRow = find(key(:,kbtrig),1);
            trig = sum(trigRow > 0);
            
            if trig
                trigTime = time(trigRow);
                key = setxor(key,key(trigRow,:),'rows','stable'); %remove trigger from keypresses
                time = setxor(time,time(trigRow));
                %trigSet = trigTime + trigOffset + fixlen; %next trial onset                
            else
                trigTime = NaN;
                %trigSet = adjBlOnset + (trial*10); %first trial onset
            end
            
            %was there a response, when, and calc rt
            
            if ~sum(key)
                respOnset = NaN;
                RT = 0;
                key = 1;
                keyPress = 0;
            else
                respOnset = time(end);
                RT = respOnset - queryOnset;
                keyPress = kb(kb==find(key(end,:),1),2);
            end
            
            
            % Record whether understood gist (ought not to)           
            
            corr_resp = 3;
            if find(key(end,:),1) == kb(3,1)
                d_acc = 1;
            else
                d_acc = 0;
            end            
            
            
            %print results     'pt,Age,Sex,CB,Block,Trial,cond,attend,speech,sent num,Sentence,mot file,targets,num Dots,exp onset,block trigger onset,block onset,trial que onset,movie onset,query onset,ITI onset,trigger,trigger time,Correct Response,key press,key time,acc,rt\n' 
            fprintf(resultsFile,'%d,%d,%s,%d,%d,%d,%d,%d,%d,%d,%s,%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%d,%d,%.4f,%d,%.4f\n', participantNumber,age,sex,cb_ans,block,-2,nan,nan,nan,nan,'rotated',motfile,numTrack,numDots,expOnset,blTrigOnset,blOnset,trOnset,movOnset,queryOnset,fixOn,trig,trigTime,corr_resp,keyPress,respOnset,d_acc,RT);
            clear locs sent            
            
            adjBlOnset = adjBlOnset + 10; %compensate timing for dummy trial
        end    

        
    end
    
    %if break block, give a break until '0' is pressed
    if ~mod(block,breakBlock) && block < nBlocks
        
        Screen('FillRect',window,backgroundColour);
        Screen('Flip',window);
        
        while ~keyCode(kbstart)
            [~,~,keyCode] = KbCheck(-1); % check for a keypress
        end
        
    end
    
end

fclose('all');

Screen('FillRect',window,backgroundColour);
Screen('TextSize',window,40);
finish = ('END');
finishW = TextBounds(window, finish);
Screen('DrawText',window,finish,xMid-(finishW(3)-finishW(1))/2,yMid, white)
Screen('Flip',window);

pause;
Screen('CloseAll');


end



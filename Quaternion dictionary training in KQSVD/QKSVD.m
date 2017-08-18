function [Dictionary,output,Dictionary1,CoefMatrix1,err] = QKSVD(...
    Data,... % an nXN matrix that contins N signals (Y), each of dimension n.
    param)
% =========================================================================
%                          K-SVD algorithm
% =========================================================================
% The K-SVD algorithm finds a dictionary for linear representation of
% signals. Given a set of signals, it searches for the best dictionary that
% can sparsely represent each signal. Detailed discussion on the algorithm
% and possible applications can be found in "The K-SVD: An Algorithm for 
% Designing of Overcomplete Dictionaries for Sparse Representation", written
% by M. Aharon, M. Elad, and A.M. Bruckstein and appeared in the IEEE Trans. 
% On Signal Processing, Vol. 54, no. 11, pp. 4311-4322, November 2006. 
% =========================================================================
% INPUT ARGUMENTS:
% Data                         an nXN matrix that contins N signals (Y), each of dimension n. 
% param                        structure that includes all required
%                                 parameters for the K-SVD execution.
%                                 Required fields are:
%    K, ...                    the number of dictionary elements to train
%    numIteration,...          number of iterations to perform.
%    errorFlag...              if =0, a fix number of coefficients is
%                                 used for representation of each signal. If so, param.L must be
%                                 specified as the number of representing atom. if =1, arbitrary number
%                                 of atoms represent each signal, until a specific representation error
%                                 is reached. If so, param.errorGoal must be specified as the allowed
%                                 error.
%    preserveDCAtom...         if =1 then the first atom in the dictionary
%                                 is set to be constant, and does not ever change. This
%                                 might be useful for working with natural
%                                 images (in this case, only param.K-1
%                                 atoms are trained).
%    (optional, see errorFlag) L,...                 % maximum coefficients to use in OMP coefficient calculations.
%    (optional, see errorFlag) errorGoal, ...        % allowed representation error in representing each signal.
%    InitializationMethod,...  mehtod to initialize the dictionary, can
%                                 be one of the following arguments: 
%                                 * 'DataElements' (initialization by the signals themselves), or: 
%                                 * 'GivenMatrix' (initialization by a given matrix param.initialDictionary).
%    (optional, see InitializationMethod) initialDictionary,...      % if the initialization method 
%                                 is 'GivenMatrix', this is the matrix that will be used.
%    (optional) TrueDictionary, ...        % if specified, in each
%                                 iteration the difference between this dictionary and the trained one
%                                 is measured and displayed.
%    displayProgress, ...      if =1 progress information is displyed. If param.errorFlag==0, 
%                                 the average repersentation error (RMSE) is displayed, while if 
%                                 param.errorFlag==1, the average number of required coefficients for 
%                                 representation of each signal is displayed.
% =========================================================================
% OUTPUT ARGUMENTS:
%  Dictionary                  The extracted dictionary of size nX(param.K).
%  output                      Struct that contains information about the current run. It may include the following fields:
%    CoefMatrix                  The final coefficients matrix (it should hold that Data equals approximately Dictionary*output.CoefMatrix.
%    ratio                       If the true dictionary was defined (in
%                                synthetic experiments), this parameter holds a vector of length
%                                param.numIteration that includes the detection ratios in each
%                                iteration).
%    totalerr                    The total representation error after each
%                                iteration (defined only if
%                                param.displayProgress=1 and
%                                param.errorFlag = 0)
%    numCoef                     A vector of length param.numIteration that
%                                include the average number of coefficients required for representation
%                                of each signal (in each iteration) (defined only if
%                                param.displayProgress=1 and
%                                param.errorFlag = 1)
% =========================================================================

if (~isfield(param,'displayProgress'))
    param.displayProgress = 0;
end
totalerr(1) = 99999;
if (isfield(param,'errorFlag')==0)
    param.errorFlag = 0;
end

if (isfield(param,'TrueDictionary'))
    displayErrorWithTrueDictionary = 1;
    ErrorBetweenDictionaries = zeros(param.numIteration+1,1);
    ratio = zeros(param.numIteration+1,1);
else
    displayErrorWithTrueDictionary = 0;
	ratio = 0;
end
if (param.preserveDCAtom>0)
    FixedDictionaryElement = 1/(2.0000001*sqrt(size(Data,1)))*ones(size(Data,1),1,4);
    %FixedDictionaryElement(:,:,1) = 0;
else
    FixedDictionaryElement = [];
end
% coefficient calculation method is OMP with fixed number of coefficients

if (size(Data,2) < param.K)
    disp('Size of data is smaller than the dictionary size. Trivial solution...');
    Dictionary = Data(:,1:size(Data,2),:);
    return;
elseif (strcmp(param.InitializationMethod,'DataElements'))
%    Dictionary(:,1:param.K-param.preserveDCAtom,:) = Data(:,2:param.K-param.preserveDCAtom+1,:);
    Dictionary(:,1:param.K-param.preserveDCAtom,:) = Data(:,1:param.K-param.preserveDCAtom,:);
elseif (strcmp(param.InitializationMethod,'GivenMatrix'))
    Dictionary(:,1:param.K-param.preserveDCAtom,:) = param.initialDictionary(1:size(Data,1),1:param.K-param.preserveDCAtom,:);
end
% reduce the components in Dictionary that are spanned by the fixed
% elements
if (param.preserveDCAtom)
    %tmpMat = FixedDictionaryElement \ Dictionary;
    %Dictionary = Dictionary - FixedDictionaryElement*tmpMat;
    for tempi = 1:param.K-1
        Dictionary(:,tempi,:) = avoidDC2(Dictionary(:,tempi,:));
    end;
else
    for tempi = 1:param.K
        len = sum(Dictionary(:,tempi,1).^2 + Dictionary(:,tempi,2).^2  + Dictionary(:,tempi,3).^2 + Dictionary(:,tempi,4).^2);
        Dictionary(:,tempi,:) = Dictionary(:,tempi,:)/sqrt(len);
    end;
end

%Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); % multiply in the sign of the first element.
%totalErr = zeros(1,param.numIteration);

% the K-SVD algorithm starts here.
err = zeros(1,param.numIteration);
for iterNum = 1:param.numIteration
    % find the coefficients
    if (param.errorFlag==0)
        %CoefMatrix = mexOMPIterative2(Data, [FixedDictionaryElement,Dictionary],param.L);
        CoefMatrix = QOMP([FixedDictionaryElement,Dictionary],Data, param.L);
    else 
        %CoefMatrix = mexOMPerrIterative(Data, [FixedDictionaryElement,Dictionary],param.errorGoal);
        CoefMatrix = QOMPerr([FixedDictionaryElement,Dictionary],Data, param.errorGoal);
        param.L = 1;
    end
    replacedVectorCounter = 0;
    rand('state',51);
	rPerm = randperm(size(Dictionary,2));
    for j = rPerm
        [betterDictionaryElement,CoefMatrix,addedNewVector] = I_findBetterDictionaryElement(Data,...
            [FixedDictionaryElement,Dictionary],j+size(FixedDictionaryElement,2),...
            CoefMatrix ,param.L);
        Dictionary(:,j,:) = betterDictionaryElement;
        if (param.preserveDCAtom)
            %tmpCoef = FixedDictionaryElement\betterDictionaryElement;
            %Dictionary(:,j) = betterDictionaryElement - FixedDictionaryElement*tmpCoef;
            %Dictionary(:,j) = Dictionary(:,j)./sqrt(Dictionary(:,j)'*Dictionary(:,j));
            Dictionary(:,j,:) = avoidDC2(Dictionary(:,j,:));
        end
        replacedVectorCounter = replacedVectorCounter+addedNewVector;
    end

    if (iterNum>1 & param.displayProgress)
        if (param.errorFlag==0)
            tempout = Qmult([FixedDictionaryElement,Dictionary],CoefMatrix);
            output.totalerr(iterNum-1) = sqrt(sum(sum(sum(Data - tempout).^2))/prod(size(Data)));
            disp(['Iteration   ',num2str(iterNum),'   Total error is: ',num2str(output.totalerr(iterNum-1))]);
        else
            output.numCoef(iterNum-1) = length(find(CoefMatrix))/size(Data,2);
            disp(['Iteration   ',num2str(iterNum),'   Average number of coefficients: ',num2str(output.numCoef(iterNum-1))]);
        end
    end
    if (displayErrorWithTrueDictionary ) 
        [ratio(iterNum+1),ErrorBetweenDictionaries(iterNum+1)] = I_findDistanseBetweenDictionaries(param.TrueDictionary,Dictionary);
        disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',num2str(ratio(iterNum+1))]));
        output.ratio = ratio;
    end
       
    
    former_Dictionary = Dictionary;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dictionary = I_clearDictionary([FixedDictionaryElement,Dictionary],CoefMatrix,Data,param.preserveDCAtom);
    Dictionary = Dictionary(:,param.preserveDCAtom+1:end,:);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(iterNum == param.numIteration)
        if(former_Dictionary == Dictionary)
            overnum = 0
        else overnum = 1
        end;
    end;
    
    if(iterNum == 1)
        Dictionary1 = [FixedDictionaryElement,former_Dictionary];
        CoefMatrix1 = CoefMatrix;
        temp_recovery1 = Qmult(Dictionary1,CoefMatrix1);
        temp_recovery1(:,:,1) = 0;
        err1 = norm(Data(:)-temp_recovery1(:),2);
    end;
    
    if (isfield(param,'waitBarHandle'))
        waitbar(iterNum/param.counterForWaitBar);
    end
    
    temp_Dictionary = [FixedDictionaryElement,former_Dictionary];
    temp_CoefMatrix = CoefMatrix;
    temp_recovery = Qmult(temp_Dictionary,temp_CoefMatrix);
    temp_recovery(:,:,1) = 0;
    err(iterNum) = norm(Data(:)-temp_recovery(:),2);
    

    if(err(iterNum) < err1)
        Dictionary1 = temp_Dictionary;
        CoefMatrix1 = temp_CoefMatrix;
        err1 = err(iterNum);
    end;
    
end

if(former_Dictionary == Dictionary)
    output.CoefMatrix = CoefMatrix;
else
    output.CoefMatrix = QOMP([FixedDictionaryElement,Dictionary],Data, param.L);
end;
Dictionary = [FixedDictionaryElement,Dictionary];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findBetterDictionaryElement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [betterDictionaryElement,CoefMatrix,NewVectorAdded] = I_findBetterDictionaryElement(Data,Dictionary,j,CoefMatrix,numCoefUsed)
if (length(who('numCoefUsed'))==0)
    numCoefUsed = 1;
end
CoefMatrix1 = CoefMatrix(:,:,1).^2 + CoefMatrix(:,:,2).^2  + CoefMatrix(:,:,3).^2 + CoefMatrix(:,:,4).^2;
relevantDataIndices = find(CoefMatrix1(j,:)); % the data indices that uses the j'th dictionary element.

if (length(relevantDataIndices)<1) %(length(relevantDataIndices)==0)
    disp('x =======================-1 ');
    ErrorMat1 = Data-Qmult(Dictionary,CoefMatrix);
    ErrorMat = ErrorMat1(:,:,1).^2 + ErrorMat1(:,:,2).^2  + ErrorMat1(:,:,3).^2 + ErrorMat1(:,:,4).^2;
    ErrorNormVec = sum(ErrorMat);
    [d,i] = max(ErrorNormVec);
    betterDictionaryElement = Data(:,i,:);%ErrorMat(:,i); 
    len = sum(betterDictionaryElement(:,:,1).^2 + betterDictionaryElement(:,:,2).^2  + betterDictionaryElement(:,:,3).^2 + betterDictionaryElement(:,:,4).^2);
    betterDictionaryElement = betterDictionaryElement./sqrt(len);
    %betterDictionaryElement = betterDictionaryElement.*sign(betterDictionaryElement(1));
    CoefMatrix(j,:,:) = 0;
    NewVectorAdded = 1;
    return;
end

NewVectorAdded = 0;
tmpCoefMatrix = CoefMatrix(:,relevantDataIndices,:); 
tmpCoefMatrix(j,:,:) = 0;% the coeffitients of the element we now improve are not relevant.
errors =(Data(:,relevantDataIndices,:) - Qmult(Dictionary,tmpCoefMatrix)); % vector of errors that we want to minimize with the new element
% % the better dictionary element and the values of beta are found using svd.
% % This is because we would like to minimize || errors - beta*element ||_F^2. 
% % that is, to approximate the matrix 'errors' with a one-rank matrix. This
% % is done using the largest singular value.

%errors(:,:,1) = 0;

[betterDictionaryElement1,singularValue,betaVector1] = qsvd(errors);
betterDictionaryElement = betterDictionaryElement1(:,1,:);
betaVector2 = betaVector1(:,1,:);
betaVector(:,:,1) = betaVector2(:,:,1)';
betaVector(:,:,2) = -betaVector2(:,:,2)';
betaVector(:,:,3) = -betaVector2(:,:,3)';
betaVector(:,:,4) = -betaVector2(:,:,4)';
CoefMatrix(j,relevantDataIndices,:) = singularValue(1,1) * betaVector;% *signOfFirstElem

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findDistanseBetweenDictionaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ratio,totalDistances] = I_findDistanseBetweenDictionaries(original,new)
% first, all the column in oiginal starts with positive values.
catchCounter = 0;
totalDistances = 0;
for i = 1:size(new,2)
    new(:,i) = sign(new(1,i))*new(:,i);
end
for i = 1:size(original,2)
    d = sign(original(1,i))*original(:,i);
    distances =sum ( (new-repmat(d,1,size(new,2))).^2);
    [minValue,index] = min(distances);
    errorOfElement = 1-abs(new(:,index)'*d);
    totalDistances = totalDistances+errorOfElement;
    catchCounter = catchCounter+(errorOfElement<0.01);
end
ratio = 100*catchCounter/size(original,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 I_clearDictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data,preserveDCAtom)
tot = 0;
similar = 0;
little_used = 0;
T2 = 0.0; % need more experience to revise
T1 = 0; % tricky choice
K=size(Dictionary,2);
Er1 = Data-Qmult(Dictionary,CoefMatrix);
Er = sum(Er1(:,:,1).^2 + Er1(:,:,2).^2 + Er1(:,:,3).^2 + Er1(:,:,4).^2);

ratio = zeros(K,K);
for x = 1:K
    for y = x+1:K
        aa = Qmult(Qpinv(Dictionary(:,x,:)),Dictionary(:,y,:));
        rem = Dictionary(:,y,:)-Qmult(Dictionary(:,x,:),aa);
        temperr = sum(rem(:,:,1).^2+rem(:,:,2).^2+rem(:,:,3).^2+rem(:,:,4).^2);
        modres = sum(Dictionary(:,y,1).^2+Dictionary(:,y,2).^2+Dictionary(:,y,3).^2+Dictionary(:,y,4).^2);
        ratio(x,y) = sqrt(temperr);
    end;
end;
ratio = ratio + ratio';
ratio = ratio + diag(ones(K,1));
CoefMatrix1 = CoefMatrix(:,:,1).^2 + CoefMatrix(:,:,2).^2  + CoefMatrix(:,:,3).^2 + CoefMatrix(:,:,4).^2;
for jj=1:1:K,
    nonzeroNo = length(find(abs(CoefMatrix1(jj,:))>1e-7));
    if min(ratio(jj,:))<T2 | nonzeroNo<=T1 ,
        tot = tot + 1;
        if min(ratio(jj,:))<T2
            similar = similar + 1;
        else little_used = little_used + 1;
        end;
        [val,pos]=max(Er);
        Er(pos(1))=0;
        if (preserveDCAtom == 1)
            Dictionary(:,jj,:) = avoidDC2(Data(:,pos(1),:));
        else
            len_D = sum(Data(:,pos(1),1).^2 + Data(:,pos(1),2).^2  + Data(:,pos(1),3).^2 + Data(:,pos(1),4).^2);
            Dictionary(:,jj,:)=Data(:,pos(1),:)/sqrt(len_D);
        end;
        ratio = zeros(K,K);
        for x = 1:K
            for y = x+1:K
                aa = Qmult(Qpinv(Dictionary(:,x,:)),Dictionary(:,y,:));
                rem = Dictionary(:,y,:)-Qmult(Dictionary(:,x,:),aa);
                temperr = sum(rem(:,:,1).^2+rem(:,:,2).^2+rem(:,:,3).^2+rem(:,:,4).^2);
                modres = sum(Dictionary(:,y,1).^2+Dictionary(:,y,2).^2+Dictionary(:,y,3).^2+Dictionary(:,y,4).^2);
                ratio(x,y) = sqrt(temperr);
            end;
        end;
        ratio = ratio + ratio';
        ratio = ratio + diag(ones(K,1));
    end;
        if jj == K
            tot
            similar
            little_used
        end;
end;

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  I_clearDictionary
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data,preserveDCAtom)
% % T2 = 0.99;
% % T1 = 3;
% % K = size(Dictionary,2);
% % temp = (Data - Qmult(Dictionary,CoefMatrix)).^2;
% % tempEr = temp(:,:,1)+temp(:,:,2)+temp(:,:,3)+temp(:,:,4);
% % Er = sum(tempEr,1);
% % DT(:,:,1)=Dictionary(:,:,1)';
% % DT(:,:,2)=-Dictionary(:,:,2)';
% % DT(:,:,3)=-Dictionary(:,:,3)';
% % DT(:,:,4)=-Dictionary(:,:,4)';
% % G = Qmult(DT,Dictionary);
% % G(:,:,1) = G(:,:,1)-diag(diag(G(:,:,1)));G(:,:,2) = G(:,:,2)-diag(diag(G(:,:,2)));G(:,:,3) = G(:,:,3)-diag(diag(G(:,:,3)));G(:,:,4) = G(:,:,4)-diag(diag(G(:,:,4)));
% % G = Q_abs(G);
% % for jj=1:1:K,
% %     if max(G(jj,:))>T2 || length(find(abs(CoefMatrix(jj,:,:))>1e-7))<=T1 ,%字典的第jj列与某一列相关性>1，或者使用字典第jj列的系数矩阵第j行的非零个数少于3
% %         [val,pos]=max(Er);
% %         Er(pos(1))=0;%避免再次用到Data的该列
% %         len_D = sum(Data(:,pos(1),1).^2 + Data(:,pos(1),2).^2  + Data(:,pos(1),3).^2 + Data(:,pos(1),4).^2);
% %         Dictionary(:,jj,:)=Data(:,pos(1),:)/sqrt(len_D);
% %         
% %         DT(:,:,1)=Dictionary(:,:,1)';
% %         DT(:,:,2)=-Dictionary(:,:,2)';
% %         DT(:,:,3)=-Dictionary(:,:,3)';
% %         DT(:,:,4)=-Dictionary(:,:,4)';
% %         G = Qmult(DT,Dictionary);
% %         G(:,:,1) = G(:,:,1)-diag(diag(G(:,:,1)));G(:,:,2) = G(:,:,2)-diag(diag(G(:,:,2)));G(:,:,3) = G(:,:,3)-diag(diag(G(:,:,3)));G(:,:,4) = G(:,:,4)-diag(diag(G(:,:,4)));
% %         G = Q_abs(G);
% %     end;
% % end;

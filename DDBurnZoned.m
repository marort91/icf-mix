function DDBurnZoned()

clf

fid=fopen('/home/mario/ICF_LANL/Vold/TeFinal.txt','w');
fid1=fopen('/home/mario/ICF_LANL/Vold/TiFinal.txt','w');
fid2=fopen('/home/mario/ICF_LANL/Vold/TradFinal.txt','w');
%fid1=fopen('/home/mario/ICF_LANL/Vold/nD.txt','w');

%Constants

JperKeV=1.6e-16; c=2.9979e8; sigma=5.67e-8; Na=6.022e23;
ccm3=1e6; KperKeV4=(1.1605e7)^4; Cv=Na*ccm3*JperKeV; Ea=820; %3520  %4850;
sigmacompton=6.653e-29; z2=1; m1=9.11e-31; e=1.602176e-19; %/1.67e-27;
kb=8.617332e-8;

%Initial Conditions

Z=10;        %Zones

%Converges for nDIC<0.65....

nDIC(1,1:Z)=1.0; nTIC(1,1:Z)=0; naIC(1,1:Z)=0;
TiIC(1,1:Z)=5; TeIC(1,1:Z)=2000; TradIC(1,1:Z)=5;
TiIC(1,:)=1.0+exp(-((1:Z)-(Z/2+0.5)).^2);

%First Time Step

nDo=nDIC; nTo=nTIC; nao=naIC;
%nDj=nDo; nTj=nTo;
ni=nDo+nTo;
ne=ni;
na=nao;

nD=nDIC;
nT=nTIC;

%ni=nDo+nTo; ne=ni; na=nao;
%nij=ni; nej=ne; naj=na;

Tio=TiIC; Teo=TeIC; Tro=TradIC;
Ti=TiIC; Te=TeIC; Trad=TradIC;

nDm=nD; nTm=nT; Tim=Ti; Tem=Te; Trm=Trad;

Ero=4.*sigma./c.*KperKeV4.*Tro.^4;
Eblackbody=4.*sigma./c.*KperKeV4.*Teo.^4;
Cvr=16.*sigma./c.*KperKeV4.*Tro.^3;

Er=Ero; Ebb=Eblackbody;

%fprintf('Hello \n ',Ero)

errormax=1e6;

N=1.0e3; iter=10; %iter=10;

time=0; tf=1e-9; dt=tf/N; tol=1e-6;

for i=1:N
    
    time=time+dt;
    error=errormax;
    
    for k=1:Z
        
        j=0;
        error=1;
        
        while (j<iter && error>tol)
            
            sigv(k)=sigvDD(Ti(k));
            rxnrateDD(k)=0.25*(nDo(k)/2)*(nDo(k)/2)*Na*ccm3*sigv(k);
            
            k1(k)=-dt*rxnrateDD(k);
            k2(k)=-dt*Na*0.25*ccm3*sigv(k)*(nDo(k)+k1(k)/2)^2;
            
            %nD=2*(nDo+k2);
            nD(k)=(nDo(k)+k2(k));
            
            fprintf(fid1,'%12.8f %12.8f %12.8f \n',nD,i,j);
            
            %nD=nDo/(1+nT*Na*ccm3*sigv*dt);
            %nT=nTo/(1+nD*Na*ccm3*sigv*dt);
            
            na(k)=(nao(k)+0.5*(nDo(k)-nD(k)));
            np(k)=na(k);
            nT(k)=(nTo(k)+0.5*(nDo(k)-nD(k)));
            
            %na=2*(nao+0.5*(nDo-nD));
            %nT=2*(nTo+0.5*(nDo-nD));
            
            ni(k)=nD(k)+nT(k)+na(k);
            ni(k)=2*ni(k);
            ne(k)=ni(k);
            
            fai(k)=0; 1/(1+32/Te(k));
            fae(k)=0; 1-fai(k);
            
            wie=4e10*ne(k)*Cv*(TeIC(k)/Te(k))^1.5;
            
            Kp=0.365*Te(k)^(-3.48447);
            rhoKp=2.5*ni(k)*Kp*100;
            
            wer=rhoKp*4*sigma*KperKeV4*(Te(k)^2+Tro(k)^2)*(Te(k)+Tro(k));
            wc=16*sigmacompton*ne*Na*ccm3*sigma*KperKeV4*Tro(k)^4/511;
            
            wer=wer+wc;
            
            Cvr=16*sigma/c*KperKeV4*Tro(k)^3;
            
            %       Ti=(Tio(k)+(Ea*fai(k)*Cv*0.5*rxnrateDD(k)+wie*Teo(k))*dt/(1.5*ni(k)*Cv))/...
            %            (1+wie*dt/(1.5*ni(k)*Cv));
            
            Tinew=(Tio(k)+(Ea*fai(k)*Cv*0.5*rxnrateDD(k)+wie*Te(k))*dt/(1.5*ni(k)*Cv))/...
                (1+wie*dt/(1.5*ni(k)*Cv));
            
            %       Te=(Teo+(Ea*fae*Cv*0.5*rxnrateDD+wie*Tio+wer*Tro)*dt/(1.5*ne*Cv))/...
            %            (1+(wie+wer)*dt/(1.5*ne*Cv));
            
            Tenew=(Teo(k)+(Ea*fae(k)*Cv*0.5*rxnrateDD(k)+wie*Ti(k)+wer*Trad(k))*dt/(1.5*ne(k)*Cv))/...
                (1+(wie+wer)*dt/(1.5*ne(k)*Cv));
            
            %       Trad=(Tro+wer*Teo*dt/Cvr)/(1+wer*dt/Cvr);
            
            Tradnew=(Tro(k)+wer*Te(k)*dt/Cvr)/(1+wer*dt/Cvr);
            
            error=max([abs(Tinew-Tim),abs(Tenew-Tem),abs(Tradnew-Trm)]);
            
          %  if error<tol
         %       break
         %   end
            
            nDm(k)=nD(k);
            nTm(k)=nT(k);
            Tim(k)=Tinew;
            Tem(k)=Tenew;
            Trm(k)=Tradnew;
            
            Tifinal(k)=Tinew;
            Tefinal(k)=Tenew;
            Trfinal(k)=Tradnew;
            
            j=j+1;
            
        end
        
        %fprintf(fid,'%12.8e %12.8f %12.8f %12.8f %12.8f %12.8f \n',time,nD,nT,Ti,Te,Trad);
        
        
        %fwrite(fid, '%12.8d %12.8d %12.8d %12.8d %12.8d %12.8d \n', time,nD,nT,Ti,Te,Trad);
        
        nDo=nD;
        nTo=nT;
        Tio=Ti;
        Teo=Te;
        Tro=Trad;
        
    end
    
    plot(Tefinal)
    pause(0.05)
    
    fprintf(fid,sprintf('%12.8e \n', Tefinal'));
    
end

fclose(fid);
fclose(fid1);

return
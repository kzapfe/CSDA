module LindbergOperators

export GaussianSmooth, DiscreteLaplacian, GaussSuavizarTemporal

#=
Module with smoothing and finite difference operators
=#

""" Gaussian function, centered """

function UnNormGauss(x,sigma)
    return exp(-x*x/(2*sigma))
end


""" Temporal Gaussian Smoothing. Second argument is sigma in
discrete frame units """

function GaussSuavizarTemporal(Datos,Sigma=3)

    medioancho=ceil(Sigma*3)
    colchon=ones(medioancho)
    result=zeros(size(Datos))
    datoscolchon=vcat(colchon*Datos[1], Datos, colchon*Datos[end])
    kernel=map(x->UnNormGauss(x,Sigma), collect(-medioancho:medioancho))
    kernel=kernel/(sum(kernel))

    #Convolution is normalized such as to preserve ratios between values.
   
    for t=medioancho+1:length(Datos)+medioancho
        result[t-medioancho]=sum(datoscolchon[t-medioancho:t+medioancho].*kernel)
    end
 
    return result
end


# Gaussian Kernel with sigma=3 electrodes
GaussianKernel=[0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
0.00002292	0.00078634	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
0.00038771	0.01330373	0.11098164	0.22508352	0.11098164	0.01330373	0.00038771
0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
0.00002292	0.00078633	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067]



""" Gausian smoothing function """
function GaussianSmooth(Datos)
    tamanodatos=size(Datos)
    result=zeros(tamanodatos)
    temp=copy(Datos)
    (mu, lu)=size(Datos)
    #problems with matrix sizes? Blame upstream developers
    arriba=reshape(temp[1,:],(1,lu))
    abajo=reshape(temp[end,:],(1,lu))
    arr3=vcat(arriba,arriba,arriba)
    aba3=vcat(abajo,abajo,abajo)
    
    temp=vcat(arr3, temp, aba3)
    
    for j=1:3
        temp=hcat(temp[:,1], temp, temp[:,end])
    end
    
    for j=4:tamanodatos[1]+3, k=4:tamanodatos[2]+3
        # Row first column second
        aux=temp[j-3:j+3,k-3:k+3]
        result[j-3,k-3]=sum(GaussianKernel.*aux)
    end
    # This convolution doesn't preserve L2 norm
    
    return result
end


#Laplace-Lindberg operator
LaplacianTerm1=[[0 1 0]; [1 -4 1]; [0 1 0]]
LaplacianTerm2=[[0.5 0 0.5]; [0 -2 0]; [0.5 0 0.5]]
LaplacianKernel=(1-1/3)*LaplacianTerm1+(1/3)*LaplacianTerm2

function DiscreteLaplacian(Datos)
    
    temp=copy(Datos)
    (mu,lu)=size(Datos)
    
    izq=reshape(temp[1,:],(1,lu))
    der=reshape(temp[end,:],(1,lu))
    # Again, padding the data
    temp=vcat(izq, temp, der)
    temp=hcat(temp[:,1], temp, temp[:,end])
    largo,ancho=size(temp)
    aux=Array{Float32}(undef, 3,3)
    result=zeros(size(temp))
    for j=2:largo-1, k=2:ancho-1
        #los indices van primero, "renglones", luego "columnas", etc
        aux=temp[j-1:j+1,k-1:k+1]
        result[j,k]=sum(LaplacianKernel.*aux)
    end
    #DO  Crop the borders
    result=result[2:end-1,2:end-1]
    return result
end


end



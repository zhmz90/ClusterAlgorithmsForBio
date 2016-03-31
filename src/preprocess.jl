

function getMutationTypes(data::Array{ASCIIString,1})
#    sort(unique(vcat(map(x->split(x,';'), data))))
    data = sort(unique(mapreduce(x->split(x,';'),vcat, data)))
    data = map(x->convert(ASCIIString,x), data)
end

function getGenes(data::Array{ASCIIString,1})
    sort(unique(data))
end

function getSampleNames(data::Array{ASCIIString,1})
    sort(unique(data))
end

function outerProduct(vector1::Array{ASCIIString,1},vector2::Array{ASCIIString,1})
    num_vector1 = length(vector1)
    num_vector2 = length(vector2)
    
    data = Array{ASCIIString,1}(num_vector1*num_vector2)
    for i = 1:num_vector1
        for j = 1:num_vector2
            data[(i-1)*num_vector2+j] = string(vector1[i],":",vector2[j])
        end
    end
    
    sort(data)
end

function preprocess(input_file)
    data = readdlm(input_file, '\t', ASCIIString)
    
    mutation_types = getMutationTypes(data[:,3])
    num_mutation_types = length(mutation_types)
    genes = getGenes(data[:,2])
    num_genes = length(genes)
    
    gene_mutations = outerProduct(genes,mutation_types)
    num_gene_mutations = length(gene_mutations)
    row_names = gene_mutations
    feature_inds = Dict{ASCIIString,Int64}(zip(gene_mutations,1:num_gene_mutations))
    
    sample_names = getSampleNames(data[:,1])
    num_samples = length(sample_names)
    col_names = sample_names
    samplename_inds = Dict{ASCIIString,Int64}(zip(sample_names, 1:num_samples))

    #TODO here need to refactor
    samplename_feature_inds = Dict{Tuple{ASCIIString,ASCIIString},Int64}()
    for sample_name in sample_names
        for gene_mutation in gene_mutations
            samplename_feature_inds[(sample_name,gene_mutation)] =
                (feature_inds[gene_mutation]-1)*(num_samples) + samplename_inds[sample_name]
        end
    end
    
    info("There are $num_samples samples")
    info("There are $num_gene_mutations features")

    nrow,ncol = size(data)
    matrix = zeros(Float64, num_samples, num_gene_mutations)
    for i = 1:nrow
        sample_name, gene_name, mutation_type = data[i,:]
        key = (sample_name,string(gene_name,":",mutation_type))

        if in(key, keys(samplename_feature_inds))
            matrix[samplename_feature_inds[key]] = 1.0
        end
    end

    
    subfeat_inds = sum(matrix,1) .!= 0
    matrix = matrix[:,subfeat_inds]
    
    writecsv("../data/matrix.csv", matrix)
    writecsv("../data/sample_names.csv", sample_names)
    @show size(sample_names)
    @show size(subfeat_inds)
    @show size(gene_mutations)
    gene_mutations = gene_mutations[subfeat_inds[:]]
    
    gene_mutations = map(gene_mutations) do gemut
        len = length(gemut)
        if len <= 20
            return gemut
        else
            return gemut[len-20:end]
        end
    end

    
    writecsv("../data/gene_mutations.csv", gene_mutations)
end

preprocess("../data/guo_gene_heatmap.txt")


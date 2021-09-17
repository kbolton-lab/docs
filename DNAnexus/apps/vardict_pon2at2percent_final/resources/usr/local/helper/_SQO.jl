import CSV
using DataFrames
import Statistics
import Combinatorics

function query(P)

    # return 0
    rename!(P,:Column1 => :CHROM)
    rename!(P,:Column2 => :POS)
    rename!(P,:Column3 => :REF)
    rename!(P,:Column4 => :ALT)
    altref = split.(P.Column6,",")
    P.altcount = [parse(Int,i[2]) for i in altref]
    P.refcount = [parse(Int,i[1]) for i in altref]
    filter!(x->x.altcount < x.refcount, P)
    P.sample_tag = ["sample" for i in P.altcount]
    filter!(x->x.altcount > 10, P)
    P.emp_vaf = P.altcount ./ (P.refcount .+ P.altcount)
    filter!(x-> length(x.REF) < 50, P)
    filter!(x-> length(x.ALT) < 50, P)
    P.Column1 = [i for i in 1:length(P[!,1])]
    return P
end

function subject(name)
    delm=","
    if name[end-2:end] == "tsv"
        delm = "\t"
    end
    M = CSV.read(name,DataFrame,header=0,delim=delm; stringtype=String)
    # return 0
    rename!(M,:Column1 => :CHROM)
    rename!(M,:Column2 => :POS)
    rename!(M,:Column3 => :REF)
    rename!(M,:Column4 => :ALT)
    altref = split.(M.Column6,",")
    M.altcount = [parse(Int,i[2]) for i in altref]
    M.refcount = [parse(Int,i[1]) for i in altref]
    M.sample_tag = ["sample" for i in M.altcount]
    M.emp_vaf = M.altcount ./ (M.refcount .+ M.altcount)
    filter!(x->x.emp_vaf > 0.001, M)
    filter!(x->x.altcount < x.refcount, M)
    M.Column1 = [i for i in 1:length(M[!,1])]
    return M
end

function part1(df,x)
    finalgroups2 = []
    perm = sortperm(x.POS) # sort by position
    x = x[perm,:] # reorder by POS
    groups = [Vector{Int32}() for i in 1:length(x.POS)]

    for row in 1:length(x.POS) # go through each variant
        temp = Vector{Int}()
        rww = row
        while true
            rww+=1
            if rww <= length(x.POS) && 0 < x.POS[rww] - x.POS[row] < 11 # see if neighbor varaints are <= 15 bp away
                push!(temp, x.Column1[rww])
            else
                break
            end
        end
        rww = row
        while true
            rww -= 1
            if rww > 0 && 0 < x.POS[row] - x.POS[rww] < 11 # see if neighbor varaints are <= 15 bp away
                push!(temp, x.Column1[rww])
            else
                break
            end
        end

        if 21 > length(temp) > 0 # if group of at least
            push!(temp, x.Column1[row])
            groups[row] = temp
        end
    end

    groups = [m for m in groups if ! isempty(m)]
    if length(groups) > 0 # see if there are any potential complex variants
        for re in 1:length(groups)
            groups[re] = sort(groups[re]) # sort each potential group by index value
        end
        if length(groups) > 1
            grps = [Vector{Int32}() for i in 1:length(groups)] # Vector{Vector{Int32}}()
            cur = groups[1]
            loc = 1
            for g in 2:length(groups) # check if there are duplicated complex
                if g == length(groups)
                    grps[loc] = groups[g]
                elseif groups[g] != cur
                    grps[loc] = cur
                    loc +=1
                    cur = groups[g]
                end
            end
            grps = [m for m in grps if ! isempty(m)]
        end
        groups = [Vector{Int32}() for i in 1:length(grps)] # Vector{Vector{Int32}}()
        tru = fill(true,length(grps))
        for g in 1:length(grps) # check if any groups are subsets of any others
            for gg in 1:length(grps)
                if tru[gg] && g != gg && issubset(grps[g], grps[gg])
                    tru[g] = false
                    @goto subsetend
                end
            end
            groups[g] = grps[g]
            @label subsetend
        end
        groups = [m for m in groups if ! isempty(m)]

        tempg = Vector{Vector{Int}}()
        # emp = df[:,:emp_vaf]
        # poscount = df[:,:altcount]
        # totalcount = df[:,:refcount]
        for g in groups
            mxx = 6 #length(g) > 6 ? 6 : length(g)
            temp = vcat([collect(Combinatorics.combinations(g,i)) for i in 1:mxx]...) #Vector{Int32}() ## gg is each starting point
            if length(temp) > 0
                append!(tempg,temp)
            end

        end
        if length(tempg) > 0
            for i in tempg
                sort!(i)
            end
            if length(tempg) > 1
                grps = [Vector{Int32}() for i in 1:length(tempg)] # Vector{Vector{Int32}}()
                cur = tempg[1]
                loc = 1
                for g in 2:length(tempg) # cut repeat small groups

                    if g == length(tempg)
                        grps[loc] = tempg[g]
                    elseif tempg[g] != cur
                        grps[loc] = cur
                        loc += 1
                        cur = tempg[g]
                    end
                end
                grps = [m for m in grps if ! isempty(m)]
            end

            groups = grps
            grps = [] #  Vector{Vector{Int32}}()
            alt = df[:,:ALT]
            pos = df[:,:POS]
            reff = df[:,:REF]
            for g in groups # group

                strt = pos[g[1]]
                complexalt = fill('-',pos[g[end]] + maximum(length.( alt[g] ) ) - strt)
                complexref = fill('-',pos[g[end]] + maximum(length.( reff[g] ) ) - strt)
                complex = fill('-',pos[g[end]] + maximum(length.( alt[g] ) ) + 30 - strt)
                #
                for gg in 1:length(g)  # element

                    for s in 1:length(alt[g[gg]]) # base
                        if complexalt[pos[g[gg]]-strt+s] == '-' || complexalt[pos[g[gg]]-strt+s] == alt[g[gg]][s]
                            complexalt[pos[g[gg]]-strt+s] = alt[g[gg]][s]
                        else
                            @goto badspot
                        end
                    end
                end
                while true
                    if complexalt[end] == '-'
                        pop!(complexalt)
                    else
                        break
                    end
                end

                for gg in 1:length(g)  # element

                    for s in 1:length(reff[g[gg]]) # base
                        if complexref[pos[g[gg]]-strt+s] == '-' || complexref[pos[g[gg]]-strt+s] == reff[g[gg]][s]
                            complexref[pos[g[gg]]-strt+s] = reff[g[gg]][s]
                        else
                            @goto badspot
                        end
                    end
                end
                while true
                    if complexref[end] == '-'
                        pop!(complexref)
                    else
                        break
                    end
                end
                complex = complexalt
                push!(grps,(g,complexalt,length(complexalt),complexref, length(complexref),complex, length(complex)  ))

                @label badspot
            end
            if length(grps) > 0
                push!(finalgroups2,((x.sample_tag[1],x.CHROM[1]),grps))
            end
        end


    end
    return finalgroups2
end

function part2(df,finalgroup,pp)

    temp=[]
    tempfinal = []

    try
        pp[( sample_tag=finalgroup[1][1], CHROM = finalgroup[1][2])]
    catch
        println(finalgroup[1][1])
        println(finalgroup[1][2])
        return []
    end
    p = pp[( sample_tag=finalgroup[1][1], CHROM = finalgroup[1][2])]

    poss = df[:,:POS]
    emp = df[:,:emp_vaf]
    reff = df[:,:REF]
    altt = df[:,:ALT]
    for j in 1:length(p[!,1])

        for k in 1:length(finalgroup[2]) # each group
            if 0 <= (p.POS[j] + length(p.ALT[j]) - poss[finalgroup[2][k][1][1]]) <= 10  # check if pindel is larger than complex and within 15 bp of start of complex
                # push!(temp, (p.Column1[j], p.REF[j], p.ALT[j], finalgroup[2][k][6], finalgroup[2][k][4], finalgroup[2][k][2], finalgroup[2][k][1], coninf[j], p.emp_vaf[j], emp[finalgroup[2][k][1]], p.CHROM[j], p.sample_tag[j], p.POS[j], poss[finalgroup[2][k][1]], reff[finalgroup[2][k][1]], altt[finalgroup[2][k][1]], p.altcount[j], string(finalgroup[2][k][6]...) ) )
                push!(temp, (p.Column1[j], p.REF[j], p.ALT[j], finalgroup[2][k][6], finalgroup[2][k][4], finalgroup[2][k][2], finalgroup[2][k][1], [], p.emp_vaf[j], emp[finalgroup[2][k][1]], p.CHROM[j], p.sample_tag[j], p.POS[j], poss[finalgroup[2][k][1]], reff[finalgroup[2][k][1]], altt[finalgroup[2][k][1]], p.altcount[j], string(finalgroup[2][k][6]...) ) )
                #               1            2           3           4                   5                        6                 7                 8            9           10                        11           12              13          14                           15                          16                          17             18                           19
            end
        end
    end
    for as in temp
        f = 0
        for i in 1:length(as[6]) # go through and make sure complex is completly within pindel and in order; in ALT
            if as[6][i] != '-'
                a=findfirst(x->x == as[6][i], as[3][f+1:end])
                if isnothing(a)
                    a=findfirst(x->x == 'N', as[3][f+1:end])
                    if isnothing(a)
                        @goto notcomplex
                    else
                        f += a
                    end
                else
                    f += a
                end
            end
        end
        f = 0
        for i in 1:length(as[5]) # go through and make sure complex is completly within pindel and in order; in REF
            if as[5][i] != '-'
                a=findfirst(x->x == as[5][i], as[2][f+1:end])
                if isnothing(a)
                    a=findfirst(x->x == 'N', as[2][f+1:end])
                    if isnothing(a)
                        @goto notcomplex
                    else
                        f += a
                    end
                else
                    f += a
                end
            end
        end
        push!(tempfinal, as)
        @label notcomplex

    end
    if length(tempfinal) == 0
        return []
    end
    compl = [i[4] for i in tempfinal]
    perm = sortperm(compl)
    return tempfinal[perm]
end


function groupsfun(df,P)
    s=groupby(df, [:sample_tag,:CHROM]) # group by sample_tag and chromosome
    finalgroups = [[] for i in 1:length(s)]
    finalgroups[1] = part1(df,s[1])
    for il in 2:length(s) # go through each group
        finalgroups[il] = part1(df,s[il])
    end
    finalgroups = vcat([m for m in finalgroups if ! isempty(m)]...)
    final = [[] for i in 1:length(finalgroups)]
    if length(final) == 0
        println("fail 2")
        return DataFrame(["nothing" "nothing";"nothing" "nothing"],:auto)
    end
    pp=groupby(P, [:sample_tag,:CHROM])

    final[1] = part2(df,finalgroups[1],pp)
    for i in 2:length(finalgroups)
        final[i] = part2(df,finalgroups[i],pp)
    end
    d=DataFrame(vcat([m for m in final if ! isempty(m)]...))
    d.literal = d[:,18]
    d.truncated = d[:,18]
    xx = d
    occurslit = zeros(Bool,length(xx[!,18]))
    occurstrunc = zeros(Bool,length(xx[!,18]))
    for i in 1:length(xx[!,18])
        xx[i,:literal] = replace(xx[i,:literal],"-"=> ".")
        if occursin(Regex(xx[i,:literal]),xx[i,3])
            occurslit[i] = true
        end
        gg = xx[i,:truncated]
        newgg = ""
        for g in 1:length(gg)
            if ! (gg[g] in ['-','.','*'])
                newgg *= gg[g]
            elseif gg[g] == '-' && ! (gg[g-1] in ['-','.','*'])
                newgg *= '.'
            elseif gg[g] == '-' && gg[g-1] == '-' && ! (gg[g-2] in ['-','.','*'])
                newgg *= '*'
            end
        end
        xx[i,:truncated] = newgg # replace(gg,"-"=> "")

        if occursin(Regex(xx[i,:truncated]),xx[i,3])
            occurstrunc[i] = true
        elseif occursin("N",xx[i,3])
            for letter in ["A","T","C","G"]
                tmp = replace(xx[i,3],"N"=>letter)
                if occursin(Regex(xx[i,:truncated]),tmp)
                    occurstrunc[i] = true
                end
            end
        end

    end
    try
        xx[1,1]
    catch
        println("fail 4")
        return DataFrame(["nothing" "nothing";"nothing" "nothing"],:auto)
    end
    dd = DataFrame(
    PindelID = xx[:,1],
    OccursTrunc = occurstrunc, OccursLiteral = occurslit,
    PindelREF = xx[:,2], PindelALT = xx[:,3],
    ComplexTrunc = xx.truncated,
    ComplexLit = xx.literal, ComplexOld = xx[:,18],
    CallerComplex = xx[:,4], CallerComplexREF= xx[:,5], CallerComplexALT= xx[:,6],
    CallerID = xx[:,7], PindelConfInt = xx[:,8], PindelVAF  = xx[:,9],
    CallerVAFmed = Statistics.median.(xx[:,10]), CallerVAF = sort.(xx[:,10]),
    CHROM  = xx[:,11], SampleTag  = xx[:,12], PindelPOS  = xx[:,13], CallerPOS  = xx[:,14],
    CallerREF= xx[:,15],CallerALT= xx[:,16] ,PindelReadSupport = xx[:,17])
    dd.ID = 1:length(dd.PindelID)
    try
        dd[1,1]
    catch
        return DataFrame(["nothing" "nothing";"nothing" "nothing"],:auto)
    end
    kep = [ ! ( (length(dd.PindelREF[i]) == 2 && length(dd.PindelALT[i]) == 1) || (length(dd.PindelREF[i]) == 1 && length(dd.PindelALT[i]) == 2) ) for i in 1:length(dd.PindelALT)]
    dd = dd[kep,:]
    select!(dd, Not([:PindelConfInt,:PindelVAF,:CallerVAF,:CallerVAFmed,:CallerID,:PindelID,:OccursTrunc,:OccursLiteral] ) )
    dd = dd[sortperm(dd.PindelReadSupport,rev=true),:]
    return unique(select(dd, Not([:ID] ) ))
end

function submain(df,p,oname,P)
    CSV.write(oname,submain2(groupsfun(df,p),P),delim="\t")
    return 0
end

function submain2(d,P)
    temp = []
    for g in groupby(d, [:PindelREF, :PindelALT, :CHROM, :SampleTag, :PindelPOS])
        temppos = []
        tempref = []
        tempalt = []
        for i in 1:nrow(g)
            push!(temppos, g[i,:CallerPOS])
            push!(tempalt, g[i,:CallerALT])
            tempalt = replace.(tempalt,"\""=>"")
            push!(tempref, g[i,:CallerREF])
            tempref = replace.(tempref,"\""=>"")
        end
        if length(temppos) < 30
            push!(temp, DataFrame(pinref=g.PindelREF[1], pinalt=g.PindelALT[1], chrom=g.CHROM[1], sample=g.SampleTag[1], pos=g.PindelPOS[1], calpos=[temppos], calref=[tempref], calalt=[tempalt]) )
        end
    end
    # P = CSV.read(qname,DataFrame,header=0; stringtype=String)
    rename!(P,:Column1 => :CHROM)
    rename!(P,:Column2 => :POS)
    rename!(P,:Column3 => :REF)
    rename!(P,:Column4 => :ALT)
    rename!(P,:Column5 => :FILTER)
    rename!(P,:Column6 => :AD)
    altref = split.(P.AD,",")
    P = leftjoin(P, select(vcat(temp...),Not(:sample) ), on=[:CHROM => :chrom, :POS => :pos, :REF => :pinref, :ALT => :pinalt] )
    P.altcount = [parse(Int,i[2]) for i in altref]
    P.refcount = [parse(Int,i[1]) for i in altref]
    P.binary = [ismissing(i) ? 0 : 1 for i in P.calpos]
    P.calpos_best = [ismissing(P[i,:calpos]) ? missing : P[i,:calpos][[ ismissing(i) ? missing : argmax(i) for i in [ ismissing(i) ? missing : [sum(length.(j)) for j in i] for i in P[:,:calpos] ]][i]] for i in 1:nrow(P)]
    P.calalt_best = [ismissing(P[i,:calalt]) ? missing : P[i,:calalt][[ ismissing(i) ? missing : argmax(i) for i in [ ismissing(i) ? missing : [sum(length.(j)) for j in i] for i in P[:,:calalt] ]][i]] for i in 1:nrow(P)]
    P.calref_best = [ismissing(P[i,:calref]) ? missing : P[i,:calref][[ ismissing(i) ? missing : argmax(i) for i in [ ismissing(i) ? missing : [sum(length.(j)) for j in i] for i in P[:,:calref] ]][i]] for i in 1:nrow(P)]
    return P
end


function main()
    delm=","
    if ARGS[2][end-2:end] == "tsv"
        delm = "\t"
    end
    P = CSV.read(ARGS[2],DataFrame,header=0,delim=delm; stringtype=String)

    submain(subject(ARGS[1]), query(copy(P)),ARGS[3], P)
    return nothing
end
main()

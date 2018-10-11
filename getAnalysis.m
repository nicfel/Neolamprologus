% get at random subset of the anophelese sequences including the x
% chromosome
clear 

system('rm -r xmls')
system('mkdir xmls')

% read in all files from the sequences folder.
loci = dir('sequences/*.fasta');

% only use 100 loci
rng(1); % reproducibility
use = randsample(length(loci),100);
loci = loci(sort(use));


loci_length = zeros(length(loci),1);

for r = 1 : 1
    f = fopen('template.xml');
    g = fopen(sprintf('xmls/neolamprologus.xml'), 'w');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line, 'insert_sequences'))
            %% print the loci
            for l = 1 : length(loci)
                fprintf(g, '\t<data id="%s" name="alignment">\n', strrep(loci(l).name, '.fa',''));
                fas = fastaread(['sequences/' loci(l).name]);
                for k = 1 : length(fas)
                    fprintf(g, '\t\t<sequence id="seq_%s_%s" taxon="%s" totalcount="4" value="%s"/>\n', strrep(loci(l).name, '.fa',''), fas(k).Header, fas(k).Header, fas(k).Sequence);
                end
                fprintf(g, '\t</data>\n');
                loci_length(l) = length(fas(1).Sequence);
            end
        elseif ~isempty(strfind(line, 'insert_init_tree'))
            %% insert starting trees
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t<tree id="Tree.t:%s" name="stateNode">\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<taxonset id="TaxonSet.%s" spec="TaxonSet">\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t<alignment idref="%s"/>\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t</taxonset>\n');
                fprintf(g, '\t\t\t\t</tree>\n');

            end

        elseif ~isempty(strfind(line, 'insert_genetreeconstrait'))
             %% insert gene tree constraints
            for l = 1 : length(loci)
%                 fprintf(g, '\t\t\t\t<distribution id="genetreeconstraint1.%s" monophyletic="true"  spec="beast.math.distributions.MRCAPrior" tree="@Tree.t:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
%                 fprintf(g, '\t\t\t\t<taxonset id="ph1.%s" spec="TaxonSet">\n', strrep(loci(l).name, '.fa',''));
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neobri" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neogra" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neomar" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neooli" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neopul" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t</taxonset>\n');
%                 fprintf(g, '\t\t\t\t</distribution>\n');
%                 fprintf(g, '\t\t\t\t<distribution id="genetreeconstraint2.%s" monophyletic="true"  spec="beast.math.distributions.MRCAPrior" tree="@Tree.t:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
%                 fprintf(g, '\t\t\t\t<taxonset id="ph2.%s" spec="TaxonSet">\n', strrep(loci(l).name, '.fa',''));
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neobri" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neogra" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neomar" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neooli" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t\t<taxon idref="neopul" spec="Taxon"/>\n');
%                 fprintf(g, '\t\t\t\t</taxonset>\n');
%                 fprintf(g, '\t\t\t\t<Uniform id="Uniform.%s" name="distr" lower="9.0" upper="20"/>\n', strrep(loci(l).name, '.fa',''));
%                 fprintf(g, '\t\t\t\t</distribution>\n');

            end               

        elseif  ~isempty(strfind(line, 'insert_kappas'))
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t<parameter id="kappa.s:%s" lower="0.0" name="stateNode">2.0</parameter>\n', strrep(loci(l).name, '.fa',''));
            end        
        elseif  ~isempty(strfind(line, 'insert_mutation_rates'))
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t<parameter id="mutationRate.s:%s" lower="0.0" name="stateNode">1.0</parameter>\n', strrep(loci(l).name, '.fa',''));
            end        
        elseif  ~isempty(strfind(line, 'insert_gene_tree'))
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t<geneTree idref="Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''));
            end 
         elseif  ~isempty(strfind(line, 'insert_gene_distribution'))
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t<distribution id="geneTree.t:%s" spec="starbeast2.GeneTreeWithMigration" ploidy="2.0" populationModel="@popModelAIM.Species" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
            end
         elseif  ~isempty(strfind(line, 'insert_kappa_prior'))
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t<prior id="KappaPrior.s:%s" name="distribution" x="@kappa.s:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<LogNormal id="LogNormalDistributionModel.kappa.%s" name="distr">\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t<parameter id="RealParameter.M.%s" estimate="false" name="M">1.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t<parameter id="RealParameter.S.%s" estimate="false" name="S">1.25</parameter>\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t</LogNormal>\n');
                fprintf(g, '\t\t\t\t</prior>\n');
            end    
         elseif  ~isempty(strfind(line, 'insert_tree_likelihood'))
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t<distribution id="treeLikelihood.%s" spec="TreeLikelihood" data="@%s" tree="@Tree.t:%s" useAmbiguities="true">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<siteModel id="SiteModel.s:%s" spec="SiteModel" mutationRate="@mutationRate.s:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t<parameter id="gammaShape.s:%s" estimate="false" name="shape">1.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t<parameter id="proportionInvariant.s:%s" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t<substModel id="hky.s:%s" spec="HKY" kappa="@kappa.s:%s">\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t\t<frequencies id="empiricalFreqs.s:%s" spec="Frequencies" data="@%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t</substModel>\n');
                fprintf(g, '\t\t\t\t\t</siteModel>\n');
                fprintf(g, '\t\t\t\t\t<branchRateModel id="StrictClock.c:%s" spec="beast.evolution.branchratemodel.StrictClockModel">\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t\t<parameter idref="strictClockRate" name="clock.rate"/>\n');
                fprintf(g, '\t\t\t\t\t</branchRateModel>\n');
                fprintf(g, '\t\t\t\t</distribution>\n');
            end           
          elseif  ~isempty(strfind(line, 'insert_down'))           
             for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t\t<down idref="Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''));
             end
          elseif  ~isempty(strfind(line, 'insert_operators'))           
             for l = 1 : length(loci)
                fprintf(g, '\t\t\t\t\t<operator id="TreeScaler.t:%s" spec="ScaleOperator" scaleFactor="0.95" tree="@Tree.t:%s" weight="3.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<operator id="TreeRootScaler.t:%s" spec="ScaleOperator" rootOnly="true" scaleFactor="0.7" tree="@Tree.t:%s" weight="3.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<operator id="UniformOperator.t:%s" spec="Uniform" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<operator id="SubtreeSlide.t:%s" spec="SubtreeSlide" size="0.002" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<operator id="Narrow.t:%s" spec="Exchange" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<operator id="Wide.t:%s" spec="Exchange" isNarrow="false" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t\t\t<operator id="WilsonBalding.t:%s" spec="WilsonBalding" tree="@Tree.t:%s" weight="15.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
             end
          elseif  ~isempty(strfind(line, 'insert_mut_scaler')) 
                fprintf(g, '\t\t\t\t<operator id="FixMeanMutationRatesOperator" spec="DeltaExchangeOperator" delta="0.75" weight="%d.0">\n', length(loci));
                for l = 1 : length(loci)
                    fprintf(g, '\t\t\t\t\t<parameter idref="mutationRate.s:%s"/>\n', strrep(loci(l).name, '.fa',''));
                end
                fprintf(g, '\t\t\t\t\t<weightvector id="weightparameter" spec="parameter.IntegerParameter" dimension="%d" estimate="false" lower="0" upper="0">', length(loci));
                for l = 1 : length(loci)
                    fprintf(g, '\t%d', loci_length(l));
                end               
                fprintf(g, '</weightvector>\n');
                fprintf(g, '\t\t\t\t</operator>\n');

        elseif  ~isempty(strfind(line, 'insert_kappa_scaler'))           
          for l = 1 : length(loci)
              fprintf(g, '\t\t\t\t<operator id="KappaScaler.s:%s" spec="ScaleOperator" parameter="@kappa.s:%s" scaleFactor="0.75" weight="1.0"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
          end                  
        elseif  ~isempty(strfind(line, 'log_tree_likelihood'))           
            for l = 1 : length(loci)
                 fprintf(g, '\t\t\t<log idref="treeLikelihood.%s"/>\n', strrep(loci(l).name, '.fa',''));
            end       
        elseif  ~isempty(strfind(line, 'log_tree_height'))           
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t<log id="TreeHeight.t:%s" spec="beast.evolution.tree.TreeHeightLogger" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
            end 
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t<log id="TreeLength.t:%s" spec="starbeast2.TreeLengthLogger" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
            end 
        elseif  ~isempty(strfind(line, 'log_kappa'))           
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t<log idref="kappa.s:%s"/>\n', strrep(loci(l).name, '.fa',''));
            end 
        elseif  ~isempty(strfind(line, 'log_mut'))           
            for l = 1 : length(loci)
                fprintf(g, '\t\t\t<log idref="mutationRate.s:%s"/>\n', strrep(loci(l).name, '.fa',''));
            end 
        elseif  ~isempty(strfind(line, 'log_gene_tree'))           
             for l = 1 : length(loci)
                fprintf(g, '\t\t<logger id="treelog.t:%s" fileName="$(tree).trees" logEvery="5000000" mode="tree">\n', strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t\t<log id="TreeWithMetaDataLogger.t:%s" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:%s"/>\n', strrep(loci(l).name, '.fa',''), strrep(loci(l).name, '.fa',''));
                fprintf(g, '\t\t</logger>\n');
             end

        else
            fprintf(g, '%s', line);
        end
    end  
    fclose('all');
end


    %% make replictes of the beast xmls
for i = 0 : 9
    f = fopen(sprintf('xmls/neolamprologus.xml'));
    g = fopen(sprintf('xmls/neolamprologus_rep%d.xml',i), 'w');
    while ~feof(f)
        line = fgets(f);
        if ~isempty(strfind(line, 'fileName="'))
            fprintf(g, strrep(line, 'fileName="', ['fileName="rep' num2str(i), '_']));
        else
            fprintf(g, line);
        end
    end
    fclose(f);
    fclose(g);
end
system('rm xmls/neolamprologus.xml');


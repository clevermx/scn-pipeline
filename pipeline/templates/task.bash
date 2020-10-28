#!/bin/bash

cd {{ AnalysisFolder }}


{% if AnalysisType == "single" and Organism == "Mus musculus" and cell_ranger and not test_mode and not panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndMus }}:/home --bind {{ WhitelistDir }}:/files --bind {{ CondaFull }}:{{ CondaFull }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif AnalysisType == "single" and Organism == "Mus musculus" and not cell_ranger and not test_mode and not panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndMus }}:/home --bind {{ CondaFull }}:{{ CondaFull }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif AnalysisType == "single" and Organism == "Homo sapiens" and cell_ranger and not test_mode and not panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndHomo }}:/home --bind {{ WhitelistDir }}:/files --bind {{ CondaFull }}:{{ CondaFull }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif AnalysisType == "single"  and Organism == "Homo sapiens" and not cell_ranger and not test_mode and not panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndHomo }}:/home --bind {{ CondaFull }}:{{ CondaFull }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif AnalysisType == "single" and Organism == "Rattus norvegicus" and cell_ranger and not test_mode and not panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndRat }}:/home --bind {{ WhitelistDir }}:/files --bind {{ CondaFull }}:{{ CondaFull }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif AnalysisType == "single" and Organism == "Rattus norvegicus" and not cell_ranger and not test_mode and not panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ IndRat }}:/home --bind {{ CondaFull }}:{{ CondaFull }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif AnalysisType == "many" and not test_mode and not bulk_like and not panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ CondaFull }}:{{ CondaFull }} --bind {{ PathToAnalysis }}:{{ PathToAnalysis }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif AnalysisType == "single" and not test_mode and not bulk_like and panglao %}
snakemake -j {{ Threads }} --use-singularity --use-conda --conda-prefix {{ PathToCondaDir }} --singularity-prefix {{ PathToSingularityDir }} --singularity-args '--bind {{ CondaFull }}:{{ CondaFull }} --bind {{ Object }}:{{ Object }} --bind {{ SummaryFile }}:{{ SummaryFile }}' --verbose
{% elif test_mode %}
snakemake -j 4 --verbose
{% endif %}
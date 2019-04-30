#!/bin/bash

#PBS -N {{ job_name }}
#PBS -P nst
#PBS -l walltime={{ time }}:00:00  
#PBS -l select={{ nodes }}:ncpus={{ ntasks }}:mem={{ mem_per_cpu }}gb
{%- if array_limit %}
#PBS -J {{array_start}}-{{ array_end }}%{{array_limit}}                                                   
{%- else %}
#PBS -J {{array_start}}-{{ array_end }}                                                                   
{%- endif %}



cd {{ execution_path }}$PBS_ARRAY_INDEX

{%- if modules_unload %}
{%- for modname in modules_unload %}
module unload {{ modname }}
{%- endfor %}
{%- endif %}
{%- if modules_load %}
{%- for modname in modules_load %}
module load {{ modname }}
{%- endfor %}
{%- endif %}

{%- if preamble %}
{{preamble}}
{%- endif %}
# Get the path to the executable; should be on user's path after the modules have been loaded.
{{ exec_path }}


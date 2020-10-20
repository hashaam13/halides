      library(dplyr)
      library(tidyr)
      library(ggplot2)
      library(data.table)
      library(scatterplot3d)
      library(tidyverse)
      
      
      args = commandArgs(trailingOnly=TRUE)
      # test if there is at least one argument: if not, return an error
      if (length(args)==0) {
        stop("At least one argument must be supplied (input file).n", call.=FALSE)
      } else if (length(args)<3) {
        # default output file
        args[3] = "out.pdf"
      }
      
      
      ###########################################
      ######## GET DATA & SET PALLETE ###########
      ###########################################
      
      df <- as_data_frame(fread(args[1]))
      df <- drop_na(df)
      cols <- c('#FAD510', '#01A08A','#FF2500', '#7A5195' )
      ###############################
      ####### Water subset###########
      ###############################
      permition = 0
      full_water <- df[df$residue_aa == 'HOH',]
      if ( nrow(full_water) != 0) {
        df_CL <- df[df$halide=='CL',]
        df_BR <-  df[df$halide=='BR',]
        df_I <-  df[df$halide=='I',]
        df_F <- df[df$halide=='F',]
        number_of_sites_CL <- length(as.vector(unlist(lapply(unique(df_CL$model_name), function(x) (table(df_CL[df_CL$model_name == x,]$id))))))
        number_of_sites_I <- length(as.vector(unlist(lapply(unique(df_I$model_name), function(x) (table(df_I[df_I$model_name == x,]$id))))))  
        number_of_sites_BR <- length(as.vector(unlist(lapply(unique(df_BR$model_name), function(x) (table(df_BR[df_BR$model_name == x,]$id)))))) 
        number_of_sites_F <- length(as.vector(unlist(lapply(unique(df_F$model_name), function(x) (table(df_F[df_F$model_name == x,]$id)))))) 
        water_CL <- df[(df$residue_aa == 'HOH') & (df$halide == 'CL'),]
        water_I <- df[(df$residue_aa == 'HOH') & (df$halide == 'I'),]
        water_BR <- df[(df$residue_aa== 'HOH') & (df$halide == 'BR'),]
        water_F <- df[(df$residue_aa== 'HOH') & (df$halide == 'F'),]
        I<- as.vector(unlist(lapply(unique(water_I$model_name), function(x) (table(water_I[water_I$model_name == x,]$id)))))
        CL <- as.vector(unlist(lapply(unique(water_CL$model_name), function(x) (table(water_CL[water_CL$model_name == x,]$id)))))
        BR <-  as.vector( unlist(lapply(unique(water_BR$model_name), function(x) (table(water_BR[water_BR$model_name == x,]$id)))))
        Fl <-  as.vector( unlist(lapply(unique(water_F$model_name), function(x) (table(water_F[water_F$model_name == x,]$id)))))
        without_water_sites_CL <- number_of_sites_CL - length(CL)
        without_water_sites_I <- number_of_sites_I - length(I)
        without_water_sites_BR <- number_of_sites_BR - length(BR)
        without_water_sites_F <- number_of_sites_F - length(Fl)
        I <- c(I,rep(0,without_water_sites_I ))
        BR <- c(BR,rep(0,without_water_sites_BR ))
        CL <- c(CL,rep(0,without_water_sites_CL ))
        Fl <- c(Fl,rep(0,without_water_sites_F ))
        water <- plyr::ldply(lapply(c('CL','BR','I'), function(x)  data.frame(water_number = eval(parse(text = x)) , halide = x )), data.frame)
        permition <-  1
      }
      ###############################
      ##### Water distribution#######
      ###############################
      if (permition == 1) {
      pdf(file = paste0(args[3], 'water distribution.pdf'))
      water %>% mutate(halide=factor(halide, levels = rev(c('I', 'BR', 'CL','Fl')))) %>% 
        ggplot(aes(water_number, fill=halide))+
        geom_histogram(col='black',  alpha=0.8, bins = 20)+
        scale_y_continuous('Count')+
        scale_fill_manual(values = cols)+
        scale_color_manual(values = cols)+
        facet_wrap(halide~., scales = 'free')+
        theme_bw()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        scale_x_continuous('Water molecules number', limits = c(-1,11))
      #dev.off()
      }
      ###############################
      ####### Water distance#########
      ###############################
      if (permition == 1) {
      pdf(file = paste0(args[3], 'water distance.pdf'))
      full_water %>% mutate(distance=as.numeric(distance)) %>% 
        mutate(halide=factor(halide, levels = rev(c('F', 'CL', 'BR', 'I')))) %>% 
        ggplot( aes(halide, distance, fill=halide))+
        geom_violin(show.legend = T, alpha=0.5, aes(col=halide))+
        geom_boxplot(width=0.3, alpha=0.7, outlier.shape = NA, show.legend = T, col='indianred4')+
        scale_fill_manual(values = rev(cols))+
        scale_color_manual(values = rev(cols))+
        scale_x_discrete('Halide')+
        scale_y_continuous('Distance')+
        coord_flip()+
        theme_bw()+
        theme(legend.position='top')
      #dev.off()
      }
      ###############################
      ######## Water angles##########
      ###############################
      if (permition == 1) {
      pdf(file = paste0(args[3], 'water angles.pdf'))
      full_water %>% mutate(angle=as.numeric(angle)) %>% 
        mutate(halide=factor(halide, levels = rev(c('F', 'CL', 'BR', 'I')))) %>% 
        ggplot( aes(halide, angle, fill=halide))+
        geom_violin(show.legend = T, alpha=0.5, aes(col=halide))+
        geom_boxplot(width=0.3, alpha=0.7, outlier.shape = NA, show.legend = T, col='indianred4')+
        scale_fill_manual(values = rev(cols))+
        scale_color_manual(values = rev(cols))+
        scale_x_discrete('Halide')+
        scale_y_continuous('Angle')+
        coord_flip()+
        theme_bw()+
        theme(legend.position='top')
        #dev.off()
      }
      
      ###########################################
      #### Coordination number distributinon ####
      ###########################################
      
      pdf(file = paste0(args[3], 'Coordination_number_distr.pdf'))
      df %>% mutate(halide=factor(halide, levels = rev(c('I', 'BR', 'CL', 'F')))) %>% 
        group_by(model_name, halide, id) %>% summarise(N = n()) %>% ungroup() %>% group_by(halide) %>% 
        ggplot(aes(N, fill=halide))+
        geom_histogram(col='black',  alpha=0.8, bins = 11)+
        scale_y_continuous('Count')+
        scale_fill_manual(values = cols)+
        scale_color_manual(values = cols)+
        facet_wrap(halide~., scales = 'free')+
        theme_bw()+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        scale_x_continuous('Coordination number', limits = c(0,11))
      dev.off()
      
      ###############################
      #### Distance distribution ####
      ###############################
      df <-  df[df$atom_type != 'WATER',]
      
      pdf(file = paste0(args[3], 'Distance_distrb.pdf'))
        df %>% mutate(halide=factor(halide, levels = rev(c('F', 'CL', 'BR', 'I')))) %>% 
          ggplot(aes(halide, distance, fill=halide))+
          geom_violin(show.legend = T, alpha=0.5, aes(col=halide))+
          geom_boxplot(width=0.3, alpha=0.7, outlier.shape = NA, show.legend = T, col='indianred4')+
          scale_fill_manual(values = rev(cols))+
          scale_color_manual(values = rev(cols))+
          scale_x_discrete('Halide')+
          scale_y_continuous('Distance')+
          coord_flip()+
          theme_bw()+
          theme(legend.position="top")
      dev.off()
      
      
      ###############################
      ###### Angle distribution #####
      ###############################
      
      pdf(file = paste0(args[3], 'Angle_distrb.pdf'))
      df %>% mutate(halide=factor(halide, levels = rev(c('F', 'CL', 'BR', 'I')))) %>% 
        filter(angle>0) %>% 
        ggplot(aes(halide, angle, fill=halide))+
        geom_violin(show.legend = F, alpha=0.5, aes(col=halide))+
        geom_boxplot(width=0.3, alpha=0.7, outlier.shape = NA, show.legend = F, col='indianred4')+
        scale_fill_manual(values = rev(cols))+
        scale_color_manual(values = rev(cols))+
        scale_x_discrete('Halide')+
        scale_y_continuous('Angle')+
        coord_flip()+
        theme_bw()
      dev.off()
      
      ###############################
      ##### fASA distribution #######
      ###############################
      
          pdf(file = paste0(args[3], 'fASA_distrb.pdf'))
          df %>% mutate(halide=factor(halide, levels = rev(c('F', 'CL', 'BR', 'I')))) %>% 
            group_by(halide, model_name, id) %>% summarise(asa=mean(asa)) %>% 
            ggplot(aes(x=halide, y=asa, fill=halide))+
            geom_violin(show.legend = F, alpha=0.5, aes(col=halide), adjust=3)+
            geom_boxplot(width=0.3, alpha=0.7, outlier.shape = NA, show.legend = F, col='indianred4')+
            scale_fill_manual(values = rev(cols))+
            scale_color_manual(values = rev(cols))+
            coord_flip()+
            theme_bw()+
            scale_x_discrete('Halide')+
            scale_y_continuous('Fractial accessible surface area')
          dev.off()
      
      ###############################
      ### Single atom frequence #####
      ###############################
      
      d <- df %>% mutate(halide=factor(halide, levels = rev(c('I', 'BR', 'CL', 'F')))) %>% group_by(halide, residue_aa, atom_name) %>%
        summarise(n = n()) %>% 
        mutate(atom = paste0(atom_name, '(', residue_aa, ')')) %>% ungroup()
      
      d_sum <- df %>% group_by(halide) %>% summarise(n_halide = n())
      d <- as_data_frame(merge(d, d_sum, by = 'halide', all.x = T))
      
      pdf(file = paste0(args[3], 'Atom_freq.pdf'))
      d %>% mutate(n_res = n/n_halide) %>% 
        group_by(halide) %>% arrange(desc(n)) %>% slice(1:10) %>%
        ggplot(aes(atom, n_res, fill=halide, col=halide))+
        geom_bar(stat = 'identity', position = 'dodge', alpha=0.8, width = 0.8)+
        facet_wrap(~halide, ncol = 4, scales = 'free')+
          scale_fill_manual(values = cols)+
          scale_color_manual(values = cols)+
        scale_y_continuous('Frequency')+
        scale_x_discrete('Atom')+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 45, hjust = 1))+
        theme(strip.background = element_blank(), strip.text.x = element_blank())+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        theme(legend.position="top")+
        coord_flip()
      dev.off()
      
      ###############################
      ### Amino acid compositions ###
      ###############################
      
      df  <- as_data_frame(fread(args[2], header = F))
      df <- drop_na(df)
      # get top 10 values
      d <- df %>% separate(V1, c('model','halide','ID', 'atoms')) %>% group_by(halide, V2) %>% summarise(n=n()) %>%
        ungroup() %>% group_by(halide) %>% arrange(desc(n)) %>% slice(1:10)
      
      # get min values for freq. occurance
      sum_d <- d %>% ungroup() %>% group_by(halide) %>% summarise(sum=sum(n), min_n=min(n)/sum)
      
      d <- as_data_frame(merge(d, sum_d))
      
      # +10% threshold for removing comp. with single matching. F was omitted.
      d_CL <- d %>% dplyr::mutate(N=n/sum) %>% filter(N>min_n*1.01, halide=='CL')
      d_BR <- d %>% dplyr::mutate(N=n/sum) %>% filter(N>min_n*1.01, halide=='BR')
      d_I <- d %>% dplyr::mutate(N=n/sum) %>% filter(N>min_n*1.01, halide=='I')
      
      pdf(file = paste0(args[3], 'Compositions.pdf'))
      bind_rows(d_CL, d_BR, d_I) %>% mutate(halide=factor(halide, levels = rev(c('F', 'CL', 'BR', 'I')))) %>% 
        ggplot(aes(x=V2, y=n/sum, fill=halide))+
        geom_bar(stat='identity', position='stack', alpha=0.8)+
        theme_bw()+
        scale_fill_manual('Halide', values = cols)+
        theme(strip.background =element_rect(fill="indianred4"))+
        theme(strip.text = element_text(colour = 'white'))+
        scale_x_discrete('Compostition')+
        scale_y_continuous('Count')+
        coord_flip()
      dev.off()
    
  
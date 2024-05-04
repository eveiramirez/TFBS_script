library(ggplot2)
library(ggpattern)
library(dplyr)
library(plyr)
library(ggpubr)
library(svglite)
library(ggrepel)

TFBS_figure<-function(index, data_list){
  subdata<-data_list[[index]]
  # Obtener los limites de la figura con respecto a las coordenadas
  gene_start<-as.numeric(subdata[which(subdata$Feature=="gene"),"Start"])
  gene_end<-as.numeric(subdata[which(subdata$Feature=="gene"),"End"])
  limits_pos<-c(gene_start, gene_end,
                as.numeric(subdata[which(subdata$Feature=="intergenic_region"),
                                   c("Start","End")]))
  
  # Al restar -1, ahora el inicio y fin no representaran nucleotidos,
  # sino el inicio y final en el eje x del feature graficado que representan
  x_min<-min(limits_pos)-1
  x_max<-max(limits_pos)
  
  # Cambiar el orden para que los TFBS aparezcan por encima y 
  # no en el fondo
  subdata_first_part<-subdata[which(subdata$Feature=="TF_binding_site"),]
  subdata_last_part<-subdata[which(subdata$Feature!="TF_binding_site"),]
  subdata<-rbind(subdata_last_part,subdata_first_part)
  
  # Restar a Start -1 para que cada posicion represente un bloque de
  # tamano 1
  # Al restar -1, ahora el inicio y fin no representaran nucleotidos,
  # sino el inicio y final en el eje x del feature graficado que representan
  subdata$Start<-subdata$Start-1
  
  # Cambiar el sentido de las coordenadas para hacer la grafica en
  # direccion 5'-3' si el gen esta en cadena reverse
  gene_strand<-subdata[which(subdata$Feature=="gene"),"Strand"]
  if(gene_strand=="-"){
    # Se intercambian los inicios y finales, y se vuelven negativos
    new_ends<-subdata$Start*(-1)
    subdata$Start<-subdata$End*(-1)
    subdata$End<-new_ends
    # Se intercambian los limites del eje x, y se vuelven negativos
    x_new_min<-x_max*(-1)
    x_max<-x_min*(-1)
    x_min<-x_new_min
  }
  
  # Obtener Tamano secuencia intergenica sin contar 5'-UTR
  five_utr_end<-subdata[which(subdata$Feature=="five_prime_UTR"),"End"]
  five_utr_start<-subdata[which(subdata$Feature=="five_prime_UTR"),"Start"]
  intergenic_reg_start<-subdata[which(subdata$Feature=="intergenic_region"),"Start"]
  int_seq_length<-five_utr_start-intergenic_reg_start

  # Se resta -1 porque anteriormente al restar -1 a todos los Start hizo que 
  # aumentara el tamano +1
  int_seq_length == int_seq_length -1
  
  # CREACION DE PUNTA DE LA FLECHA
  # Obtener inicio y fin de las coordenadas del eje x del gen
  gene_start<-as.numeric(subdata[which(subdata$Feature=="gene"),"Start"])
  gene_end<-as.numeric(subdata[which(subdata$Feature=="gene"),"End"])
  # Obtener tamano del gen
  gene_length<-gene_end-gene_start
  # Obtener el tamano de 1/8 del gen, y eliminar decimales si se obtuvieron
  arrow_head_length<-trunc(gene_length/8)
  # Determinar la posiciones finales e iniciales de la punta de la flecha
  # en el eje x
  arrow_head_xmin<-gene_end-arrow_head_length
  arrow_head_xmax<-gene_end
  
  # Obtener inicio y fin de las coordenadas del eje x del 3'-UTR
  three_utr_start<-subdata$Start[which(subdata$Feature=="three_prime_UTR")]
  three_utr_end<-subdata$End[which(subdata$Feature=="three_prime_UTR")]
  # Cambiar final del gen
  subdata[which(subdata$Feature=="gene"),"End"]<-arrow_head_xmin
  
  # Generar coordenadas para la punta de la cabeza
  # Evaluar si el inicio del 3'-UTR esta despues del inicio de la punta
  if(three_utr_start>arrow_head_xmin){
    len_min<-three_utr_start-arrow_head_xmin
    len_max<-arrow_head_xmax-arrow_head_xmin
    y_min<-(len_min/len_max)*0.5
    y_max<-1-y_min
    x_positions<-c(arrow_head_xmin,arrow_head_xmin,three_utr_start,
                   three_utr_start,three_utr_start,three_utr_start,
                   arrow_head_xmax)
    y_positions<-c(0,1,y_min,y_max,y_min,y_max,0.5)
    arrowhead_positions <- data.frame(
      id = c(rep("gene",4), 
             rep("three_prime_UTR",3)),
      x_pos = x_positions,
      y_pos = y_positions
    )
    chulls <- ddply(arrowhead_positions, .(id), 
                    function(arrowhead_positions){
                      arrowhead_positions[chull(arrowhead_positions$x_pos, 
                                                arrowhead_positions$y_pos),]})
    arrowhead_positions<-chulls
    subdata<-subdata[which(subdata$Feature!="three_prime_UTR"),]
  } else{
    x_positions<-c(arrow_head_xmin,arrow_head_xmin,arrow_head_xmax)
    y_positions<-c(0,1,0.5)
    arrowhead_positions <- data.frame(
      id = "three_prime_UTR",
      x_pos = x_positions,
      y_pos = y_positions
    )
    if(three_utr_start==arrow_head_xmin){
      subdata<-subdata[which(subdata$Feature!="three_prime_UTR"),]
    } else{
      subdata[which(subdata$Feature=="three_prime_UTR"),"End"]<-arrow_head_xmin
      subdata[which(subdata$Feature=="gene"),"End"]<-three_utr_start
    }
  }
  
  # Cambiar el inicio del gen, para representar la region codificante 
  subdata[which(subdata$Feature=="gene"),"Start"]<-five_utr_end
  
  # Cambiar el nombre del feature de los TFBS de acuerdo a su cadena
  subdata$Feature[which(subdata$Strand=="+" & subdata$Feature=="TF_binding_site")]<-"TF_binding_site_plus"
  subdata$Feature[which(subdata$Strand=="-" & subdata$Feature=="TF_binding_site")]<-"TF_binding_site_m"
  
  # Cambiar la altura de la figura de acuerdo al index de la figura
  arrowhead_positions$y_pos<-arrowhead_positions$y_pos+2*(index-1)
  subdata$y_min<-subdata$y_min+2*(index-1)
  subdata$y_max<-subdata$y_max+2*(index-1)
  
  # Restar el inicio de reg intergenica para que todos los features
  # inicien en 0
  arrowhead_positions$x_pos<-arrowhead_positions$x_pos-x_min
  subdata$Start<-subdata$Start-x_min
  subdata$End<-subdata$End-x_min
  x_max<-x_max-x_min
  x_min<-0
  
  # Obtener el nuevo inicio del 5'-UTR
  five_utr_start<-subdata[which(subdata$Feature=="five_prime_UTR"),"Start"]
  # Obtener los inicios y finales
  tfbss_indices<-subdata$Feature %in% c("TF_binding_site_plus",
                                          "TF_binding_site_m")
  
  tfbss_distance<-subdata[tfbss_indices, c("Start","End")]
  
  # Obtener distancia de TFBSs del inicio del 5'-UTR
  # Si la distancia es cero, o si el inicio del 5'-UTR esta contenido en el TFBS,
  # la distancia permanece como cero
  tfbss_distance$Distance<-0
  tfbss_distance$x_cord<-five_utr_start
  tfbss_distance_xcord<-sapply(1:nrow(tfbss_distance), function(x){
    if(tfbss_distance$End[x]<five_utr_start){
      print
      return(c(-(five_utr_start-tfbss_distance$End[x]),tfbss_distance$End[x]))
    } else if(tfbss_distance$Start[x]>five_utr_start){
      return(c(tfbss_distance$Start[x]-five_utr_start,tfbss_distance$Start[x]))
    } else(return(c(0,five_utr_start)))
  })
  
  tfbss_distance$Distance<-tfbss_distance_xcord[1,]
  tfbss_distance$x_cord<-tfbss_distance_xcord[2,]
  tfbss_distance$y_cord<-1+2*(index-1)

  # Devolver los nuevos datos y tamano de reg intergenica
  subdata_list<-list(subdata,arrowhead_positions,int_seq_length,tfbss_distance)
}

# FUNCION QUE GENERA FIGURA DE LOS FEATURES DE UN ARCHIVO GFF
get_TFBS_figures<-function(directory, y_labels=NULL,title="",text_size=11){
  # Leer el archivo GFF
  inputdata<-read.table(directory, sep="\t", header=F)
  # Agregar los nombres de las columnas
  colnames(inputdata)<-c("Feature", "Start", "End", "Strand")
  
  # Extraer unicamente los features de los genes, UTRs, TFBSs, y regiones
  # intergenicas upstream
  # El programa no determina si es upstream o downstream, solo que exista
  # un feature que sea region intergenica
  features_indices<-inputdata$Feature %in% c("five_prime_UTR","three_prime_UTR",
                                             "gene","TF_binding_site",
                                             "intergenic_region")
  data<-inputdata[features_indices,]

  # Calcular la altura de los bloques dada la cadena en la que se
  # encuentre el feature
  data$y_min<-0
  data$y_max<-1
  
  # Volver los inicios y finales numericos
  data$Start<-as.numeric(data$Start)
  data$End<-as.numeric(data$End)
  
  # Obtener las posiciones de las regiones intergenicas en la tabla
  int_index<-c(which(data$Feature=="intergenic_region"), nrow(data))
  
  # Dividir los datos de acuerdo a las posiciones de las regiones intergenicas
  # Se asume que el primer feature siempre es una region intergenica, y que
  # los demas features que correspondan al gen o region intergenica estan antes
  # de la siguiente region intergenica del siguiente gen
  ind_end<-length(int_index)-1
  if(ind_end==1){
    data_list<-list(data)
  } else{
    data_list<-lapply(1:ind_end, function(x){
      if(x<ind_end){
        return(data[int_index[x]:(int_index[x+1]-1),])
      } else{
        return(data[int_index[x]:(int_index[x+1]),])
      }
    })
  }
  # Generar indices de cada elemento de la lista
  index<-1:length(data_list)
  # Obtener los datos para las figuras de cada set de datos en una lista
  figures_list<-lapply(index, TFBS_figure, data_list=data_list)
  
  # Obtener los datos de los features de cada una de las partes
  features<-lapply(index, function(x, figures){
    return(figures[[x]][[1]])
  },figures=figures_list)
  # Unir los datos
  features<-bind_rows(features, .id = "column_label")
  
  # Obtener los datos de las puntas de las flechas de los genes
  arrow_heads<-lapply(index, function(x, figures){
    return(figures[[x]][[2]])
  },figures=figures_list)
  
  # Obtener los tamanos de las regiones intergenicas
  intergenic_lengths<-sapply(index, function(x, figures){
    return(figures[[x]][[3]])
  },figures=figures_list)
  
  # Obtener las distancias del los TFBSs del inicio del 5'-UTR
  tfbss_distances<-lapply(index, function(x, figures){
    return(figures[[x]][[4]])
  },figures=figures_list)
  
  # Obtener coordenadas del extremo de cada gen
  arrow_heads_xmax<-sapply(index, function(x, arrow_heads){
    return(max(arrow_heads[[x]]$x_pos))
  },arrow_heads=arrow_heads)
  
  # Obtener extremos de la figura
  x_min<-min(features$Start)
  x_max<-max(arrow_heads_xmax)
  y_min<-min(features$y_min)
  y_max<-max(features$y_max)
  
  # Generar breaks del eje y
  y_breaks<-sapply(index, function(x){2*(x-1)})
  # Generar labels del eje y si no existen
  if(is.null(y_labels)){
    y_labels<-index
  }
    
  # Generar paleta de colores y aisgnar colores
  # Se genera la paleta de colores para los features
  features_palette<-c("#EEDD88","#99DDFF","#DDDDDD","#77AADD","#AAAA00",
                      "#FFAABB")
  
  # Generar y guardar el plot
  final_plot<-ggplot() +
    # Generar las cajas de cada feature
    geom_rect(data=features, aes(xmin=Start,xmax=End,ymin=y_min,ymax=y_max,
                                 fill=Feature, colour=Feature)) +
    # Establecer las posiciones a graficar del eje x
    # Anadir los limites -1 para que los ticks de los extremos aparezcan
    scale_x_continuous(limits=c(x_min-1, x_max+1),
                       expand = expansion(mult = c(0,.05))) +
    # Establecer los limites del eje y
    # NOTA: PARA ANADIR EL ZOOM A LA FIGURA, SE RESTO -1 A y_min
    scale_y_continuous(limits=c(y_min-0.2-1, y_max+0.2), 
                       breaks=0.5+y_breaks, 
                       labels=y_labels,
                       expand = expansion(mult = c(0.1,0.1))) +
    # Asignar colores de relleno y contorno, y las etiquetas de estos
    scale_fill_manual(values=features_palette,
                      labels=c("5'-UTR","ORF","Intergenic region",
                               "G2-TFBS (rev strand)",
                               "G2-TFBS (fwd strand)","3'-UTR")) +
    scale_colour_manual(values=features_palette) +
    # Eliminar la leyenda del color
    guides(colour="none") +
    theme(axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(),
          panel.background=element_rect(fill="white"),
          axis.text.x=element_blank(),
          axis.text.y=element_text(family="sans",size=text_size,colour="red",
                                   hjust = 0),
          legend.text=element_text(family="sans",size=text_size,colour="red"),
          legend.title=element_text(family="sans",size=text_size,colour="purple"),
          plot.title = element_text(family="sans",size=text_size,colour="red")) +
    # Generar los labels
    labs(x="",y="", fill="Features", title=title)
  
  # Anadir las cabezas de las flechas de forma separada, porque si no se hace
  # asi se tienen que separar en grupos, lo cual hara que sea necesario 
  # asignar colores a cada grupo
  for (i in index) {
    final_plot<-final_plot+
      geom_polygon(data=arrow_heads[[i]],
                   aes(x=x_pos, y=y_pos, fill=id, colour=id, group=id))+
      #Generar una linea representando el eje x
      geom_text(data=data.frame(x_cord=x_min,y_cord=0+2*(i-1)),
                aes(x=x_cord,y=y_cord),
                label=paste(round(intergenic_lengths[i]/1000,digits=1),"kb"),
                hjust="left",
                nudge_y=-0.16,family="sans", colour="blue", 
                size=text_size/.pt)+
      geom_text(data=tfbss_distances[[i]],
                aes(x=x_cord,y=y_cord,label=Distance),
                hjust="center",
                nudge_y=0.16,family="sans", colour="blue", 
                size=text_size/.pt,
                check_overlap=TRUE)
  }
  # Devolver plot
  final_plot
}


# GENERAR IMAGEN DE PRUEBA
directory<-c("input/test_dataset.tsv")
figures<-get_TFBS_figures(directory)
figures+geom_vline(xintercept=0:60)

# GENERAR IMAGEN DE GENES DE M.polymorpha SIN ZOOM DE TFBSs
directory<-c("input/Mpolymorpha_6.1.tsv")

y_labels<-c("GCAM1","WIP","NOP1")
title<-""
figures<-get_TFBS_figures(directory, y_labels, title,text_size=12)
figures


# GENERAR ZOOM DE TFBSs SOBRELAPADOS EN GCAM1
zoom_tfbs<-data.frame(Start=c(3828,5032,5333),
                      End=c(9547,8042,8343),
                      y_min=-1,
                      y_max=-0.25,
                      color=c("int","-","+"))
zoom_tfbs_pattern<-data.frame(Start=5333,
                      End=8042,
                      y_min=-1,
                      y_max=-0.25,
                      color=c("-"))

# NOTA: Si el codigo es ejecutado en RStudio, tener en cuenta que pueden
# ocurrir errores si se desea visualizar las figuras por tener un tamano
# pequeno de pantalla para el visualizar de los Plots
figure<-figures+geom_line(data=data.frame(x=c(6678,3828,6697,9547),y=c(0,-0.25,0,-0.25),
                             group=c(1,1,2,2)),aes(x=x,y=y,group=group)) +
  geom_rect(data=zoom_tfbs, aes(xmin=Start,xmax=End,ymin=y_min,ymax=y_max),
            colour=c("#DDDDDD","#77AADD","#AAAA00"),
            fill=c("#DDDDDD","#77AADD","#AAAA00")) +
  geom_rect_pattern(data=zoom_tfbs_pattern,
                    aes(xmin=Start,xmax=End,ymin=y_min,ymax=y_max),
                    pattern="stripe",
                    pattern_colour="#77AADD",
                    pattern_fill="#77AADD",
                    pattern_frequency=0.5,
                    colour="#77AADD",
                    alpha=0)
  
figure



# GENERAR IMAGEN DE GENES DE M.polymorpha SIN ZOOM DE TFBSs, del articulo
# de copas de 2023
directory<-c("input/copas_scan_result.tsv")

y_labels<-c("MAX2","SMXL")
title<-""
figure<-get_TFBS_figures(directory, y_labels, title,text_size=12)
figure


# GUARDAR IMAGEN DEL SET DE PRUEBAS
directory<-c("figuras")
svg(directory, width=12,height=1.5*10,pointsize=12)
figures
dev.off()

# GUARDAR IMAGEN DEL SET DE GENES
ggsave(filename="figures_final_v1.svg", plot=figure,
       path=directory,
       width = 30.48, height=11.43, units ="cm", dpi=1200)

# GUARDAR IMAGEN DEL SET DE GENES de COPAS
ggsave(filename="figures_copas.svg", plot=figure,
       path=directory,
       width = 30.48, height=11.43, units ="cm", dpi=1200)

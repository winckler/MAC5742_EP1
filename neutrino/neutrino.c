/* EP1 - MAC5742 - 2011/2 - #USP 3313359      */
/* Copyright (C) 2011 Gabriel A. von Winckler */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* Estrutura de pastilha individual
   mantem a vizinhanca e os calculos intermediarios
   OBS: a temperatura ee mantida em um (dois) vetor separado */
typedef struct {
    double dl;            /* altura efetiva */
    int bursted;          /* se esta estourada ou nao */
    
    double delta_atrito;  /* valor intermediario do atrito */
    double delta_diss;    /* valor intermediario da dissipacao */
    
    double *r_tile_temp;  /* temperatura pastilha a direita */
    double *l_tile_temp;  /* temperatura pastilha a esquerda */
    double *t_ring_temp;  /* temperatura anel superior */
    double *b_ring_temp;  /* temperatura anel inferior */
} tile_t;

/* Estrutura de anel
   mantem a referencia a primeira e ultima pastilha
   OBS: a temperatura ee mantida em vetor separado */
typedef struct {
    unsigned long int first_tile;
    unsigned long int last_tile;
} ring_t;

/* funcao que formata a impressao */
double print_temp(double temp, unsigned char bursted) {
    if (bursted)
        return -temp;
    else
        return temp;
}

/* funcao que executa a rotacao do sistema de coordenadas
   solitario ao elipsoide para a coordenada global.
   Essa implementacao esta de acordo com a feita na implementacao de referencia */
void rotate_system(double *p_x, double *p_y, double *p_z, 
                   double *v_x, double *v_y, double *v_z) {

    /* Rotate systems */
    double mod_p = sqrt(pow(*p_x, 2) + pow(*p_y, 2) + pow(*p_z, 2));
    *p_x /= mod_p;
    *p_y /= mod_p;
    *p_z /= mod_p;

    double a, b;
    if (*p_x == 0)
        a = M_PI / 2.0;
    else
        a = -1 * atan( *p_y / *p_x ); 
    
    double xl, yl, zl;
    xl = *v_x * cos(a) - *v_y * sin(a);
    yl = *v_x * sin(a) + *v_y * cos(a);
    
    *v_x = xl;
    *v_y = yl;
    
    xl = *p_x * cos(a) - *p_y * sin(a);
    yl = *p_x * sin(a) + *p_y * cos(a);
    
    *p_x = xl;
    *p_y = yl;  

    if (*p_z == 0)
        b = M_PI / 2.0;
    else
        b = -1 * atan( *p_x / *p_z ); 
        
    xl = *v_z * sin(b) + *v_x * cos(b);
    zl = *v_z * cos(b) - *v_x * sin(b);
    
    *v_x = xl;
    *v_z = zl;
    
    xl = *p_z * sin(b) + *p_x * cos(b);
    zl = *p_z * cos(b) - *p_x * sin(b);
    
    *p_x = xl;
    *p_z = zl;
    
    if (*p_z > 0) {
        *p_x = -*p_x;
        *p_y = -*p_y;
        *p_z = -*p_z;
    
        *v_x = -*v_x;
        *v_y = -*v_y;
        *v_z = -*v_z;
    }
}

int main(void)
{
    
    /* parametros de entrada */
    double h, a, d, alpha, t_0, delta, theta_0, theta_crit;
    double p_x, p_y, p_z;
    double v_x, v_y, v_z;
    unsigned long int steps;

    unsigned long int t, s, r; /* loop counters */
    double l; /* passo da altura */

    /* leitura do stdin */
    int read = 0;
    read += scanf("%lf%lf", &h, &a);
    read += scanf("%lf", &d);
    read += scanf("%lf%lf", &alpha, &t_0);
    read += scanf("%lf", &delta);
    read += scanf("%lf%lf", &theta_crit, &theta_0);
    read += scanf("%lf%lf%lf", &p_x, &p_y, &p_z);
    read += scanf("%lf%lf%lf", &v_x, &v_y, &v_z);
    read += scanf("%lu", &steps);

    /* verifica se tem todos os valores */
    if (read != 15) {
        printf("Error in input\n");
        return 1;
    };

    /* determina a altura inicial e o passo */
    double L = a * pow( 3 * d / M_PI , 2);
   
    /* calcula numero de aneis */
    unsigned long int n_rings = ((int) floor(h / L) - 1);
    
    /* calcula numero de pastilhas */
    unsigned long int n_tiles = 0;
    for (l = L; (l+L) < h; l += L) {
        n_tiles += (int) (2 * M_PI * sqrt(l / a)) / d;
    };
    
    printf("n_rings = %lu\n", n_rings);
    printf("n_tiles = %lu\n", n_tiles);

    /* lista de aneis e pastilhas */
    ring_t *rings = (ring_t *) malloc(n_rings * sizeof(ring_t));
    tile_t *tiles = (tile_t *) malloc(n_tiles * sizeof(tile_t));

    /* vetor com a temperatura dos aneis */
    double *rings_temp = (double *) malloc((n_rings+1) * sizeof(double));
    
    /* vetor de temperaturas (passo anterior e proximo)
       OBS: fazer fora da struct permite atulizar com um unico memcpy */
    double *tiles_cur_temp = (double *) malloc((n_tiles+1) * sizeof(double));
    double *tiles_next_temp = (double *) malloc((n_tiles+1) * sizeof(double));

    /* verifica de uma vez se conseguiu alocacao de memoria. nao ee bonito, mas... */
    if ((rings == NULL) ||
          (tiles == NULL) ||
          (rings_temp == NULL) || 
          (tiles_cur_temp == NULL) ||
          (tiles_next_temp == NULL)) {
        printf("Memory allocation error\n");
        return 2;   
    }    
    
    /* executa a rotacao no sistema de coordenadas */
    rotate_system(&p_x, &p_y, &p_z,
                  &v_x, &v_y, &v_z);

    /* inicializa a calota */
    double cover_temp = theta_0;
    double cover_vn = v_x * p_x + v_y * p_y + v_z * p_z;
    
    double cover_delta_atrito = 0;
    if (cover_vn > 0)
        cover_delta_atrito = alpha * cover_vn;

    double cover_delta_diss = delta * fabs(cover_vn);
        
    int cover_bursted = 0;

  
    unsigned long int ring = 0; /* contadores */  
    unsigned long int tile = 0;

    /* loop que gera os aneis */    
    for (l = L; (l+L) < h; l += L) {
        /* pastilhas desse anel */
        int n = (int) (2 * M_PI * sqrt(l / a)) / d;

        /* salva os limites */
        rings[ring].first_tile = tile;
        rings[ring].last_tile = tile + n - 1; 

        /* cria a geomeria */
        double z0 = l;
        double z1 = l+L;
            
        double r0 = sqrt(z0 / a);
        double r1 = sqrt(z1 / a);
        
        double alpha0 = 2 * asin( (d/2) / r0 );
        double alpha1 = 2 * asin( (d/2) / r1 );
            
        double beta0 = ((2 * M_PI) - (n * alpha0)) / n;
        double beta1 = ((2 * M_PI) - (n * alpha1)) / n;

        /* anel anterior */
        double *b_ring;
        if (ring == 0)
            b_ring = &cover_temp;
        else
            b_ring = &rings_temp[ring-1];

        /* anel posterior */
        double *t_ring;
        if (ring == n_rings - 1)
            t_ring = &rings_temp[ring];
        else
            t_ring = &rings_temp[ring+1];

        /* loop que cria cada pastilha */
        for (t=0; t < n; t++) {
            /* geometria, e mais gometria */
            double x0 = r0 * cos( t * (alpha0 + beta0) );
            double y0 = r0 * sin( t * (alpha0 + beta0) );

            double X0 = r0 * cos( (t+1) * (alpha0 + beta0) );
            double Y0 = r0 * sin( (t+1) * (alpha0 + beta0) );

            double x1 = r1 * cos( t * (alpha1 + beta1) );
            double y1 = r1 * sin( t * (alpha1 + beta1) );

            double X1 = r1 * cos( (t+1) * (alpha1 + beta1) );
            double Y1 = r1 * sin( (t+1) * (alpha1 + beta1) );
            
            /* temperatura inicial */
            tiles_cur_temp[tile] = theta_0;
            tiles[tile].dl = sqrt( pow(x0 - x1, 2) +
                                   pow(y0 - y1, 2) + 
                                   pow(z0 - z1, 2) );
            
            double p0 = X0 - x0;
            double p1 = Y0 - y0;
            double p2 = z0 - z0;
            
            double q0 = X1 - x0;
            double q1 = Y1 - y0;
            double q2 = z1 - z0;
            
            double n_x = (p1 * q2) - (p2 * q1);
            double n_y = (p2 * q0) - (p0 * q2);
            double n_z = (p0 * q1) - (p1 * q0);

            double mod_n = sqrt(pow(n_x, 2) + pow(n_y, 2) + pow(n_z, 2));
 
            double vn = (v_x * n_x + v_y * n_y + v_z * n_z) / mod_n;

            /* ja salva o valor preprocessado */
            if (vn > 0)
                tiles[tile].delta_atrito = alpha * vn;
            else
                tiles[tile].delta_atrito = 0;
            tiles[tile].delta_diss = delta * fabs(vn);
            tiles[tile].bursted = 0;

            /* Vizinhos */
            tiles[tile].r_tile_temp = &tiles_cur_temp[tile+1];
            tiles[tile].l_tile_temp = &tiles_cur_temp[tile-1];
            
            tiles[tile].t_ring_temp = t_ring;
            tiles[tile].b_ring_temp = b_ring;
            
            tile++;
        }

        /* corrigir a vizinhança das pastilhas (fecha o circulo) */
        tiles[rings[ring].first_tile].l_tile_temp = &tiles_cur_temp[rings[ring].last_tile];
        tiles[rings[ring].last_tile].r_tile_temp = &tiles_cur_temp[rings[ring].first_tile];

        /* temperatura inicial */
        rings_temp[ring] = theta_0;
        ring++;
    };    

    
    /* Main Loop: vamos ao que interessa */
    /* para cada passo... */
    for(s=0; s < steps; s++) {
        /* OBS: O calculo mais caro da formula ee o atan, mas que independe
           da pastilha, so do tempo. Por isso estou calculando fora do loop e salvando */
    	double atan_atrito = atan(pow((s+1) - t_0, 2));

        /* principal loop para paralelizar. todas as pastilhas podem ser feitas juntas */
    	#pragma omp parallel for default(shared) schedule(static, n_tiles/omp_get_num_procs())
        for(t=0; t < n_tiles; t++) {
            double perm_temp, new_temp;
         
            /* temperatura na vizinhança */
            perm_temp = ((*tiles[t].l_tile_temp + *tiles[t].r_tile_temp) * tiles[t].dl +
                         (*tiles[t].t_ring_temp + *tiles[t].b_ring_temp) * d             ) / 
                         (2 * (tiles[t].dl + d));
            
            new_temp = perm_temp;
            /* verifica se a pastilha ja estava estourada */
            if (! tiles[t].bursted) {
                new_temp = perm_temp + (tiles[t].delta_atrito * atan_atrito) - tiles[t].delta_diss;
                
                /* verifica se acabou de estourar */
                if (new_temp > theta_crit) {
                    new_temp = perm_temp;
                    tiles[t].bursted = 1;
                }
            }
            
            /* salva no vetor */
            tiles_next_temp[t] = new_temp;
        };

        /* Calcula a calota, procedimento analogo as pastilhas */
        cover_temp = rings_temp[0];
        if (! cover_bursted) {
                cover_temp =  rings_temp[0] + (cover_delta_atrito * atan_atrito) - cover_delta_diss;
            
                if (cover_temp > theta_crit) {
                    cover_temp = rings_temp[0];
                    cover_bursted = 1;
                }
        };
        
        /* Atualiza as temperaturas - so um memcpy */
        memcpy(tiles_cur_temp, tiles_next_temp, n_tiles * sizeof(double));
    
        double total = 0;
        /* Faz as médias dos aneis (e a media total para impressao) */
       	#pragma omp parallel for default(shared) reduction(+:total) schedule(dynamic, n_rings / omp_get_num_procs() / 4)
        for (r=0; r < n_rings; r++) {
            double sum = 0;
            for(t=rings[r].first_tile; t <= rings[r].last_tile; t++)
                sum += tiles_cur_temp[t];
                
            rings_temp[r] = sum / (rings[r].last_tile - rings[r].first_tile + 1);
            total += sum;
        };
        
        printf("step = %lu, average = %.6lf\n", s, total / n_tiles);
    
    };
    /* ja acabou... */
 
 
    /* arquivo de saida */
    FILE *out = fopen("neutrino.out", "w");
    fprintf(out, "%.6lf %.6lf %.6lf\n", a, h, d);
    fprintf(out, "%.6lf\n", print_temp(cover_temp, cover_bursted));
    fprintf(out, "%lu\n", n_rings);
    for (r=0; r < n_rings; r++) {
        fprintf(out, "%lu ", rings[r].last_tile - rings[r].first_tile + 1);
        for(t=rings[r].first_tile; t <= rings[r].last_tile; t++)
            fprintf(out, "%.6lf ", print_temp(tiles_cur_temp[t], tiles[t].bursted));
        fprintf(out, "\n");
    };
    fclose(out);
    
    printf("Done\n");
    return 0;
}

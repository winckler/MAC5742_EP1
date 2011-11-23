#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Esboco do EP de MAC5742 - 2/2011
# 
# Copyright (C) 2011 Gabriel A. von Winckler
#                    Max Rosan
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# version 0.3
#
# Bugs:
#  - A os pontos e a normal está calculada nas coordenadas do elipsoide

import sys
import re
import math

# Carrega o módulo 3D caso o módulo vtk estaja disponível
# no ubuntu: # apt-get install python-vtk
try:
    import vtk
except ImportError:
    print "You don't have vtk module. No visualization for you."
    vtk = False

class tile():
    # Classe das pastilhas
    def __init__(self, inputs, a, b, c, d, parent):

        # Tempo inicial e contador para incrementar (necessário para a equação do atrito)
        self.t_0 = inputs["t_0"]
        self.t = self.t_0

        # Parâmetros de entrada para as equações de atrito e dissipação
        self.alpha = inputs["alpha"]
        self.delta = inputs["delta"]

        # Temperatura inicial e temperatura em que a pastilha explode
        self.last_temp = inputs["theta_0"]
        self.max_temp = inputs["theta_crit"]

        # Vetores da velocidade e posição
        self.vel = inputs["vel"]
        self.pos = inputs["pos"]

        # Altura e base da pastilha
        self.dl = math.sqrt( (a[0]-c[0]) ** 2 + (a[1]-c[1]) ** 2 +  (a[2]-c[2]) ** 2 )
        self.d = inputs["d"]

        # Estado da pastilha
        self.bursted = False # True quando a pastilha estourou quando estourou

        # Salva os vértices para a exibição
        self.edges = (a, b, c, d)

        # produto vetorial n = (p2 - p1) x (p3 - p1)
        p = [x - y for x, y in zip (b, a)]
        q = [x - y for x, y in zip (d, a)]

        # normal na coordenada do paraboloide.
        self.normal = tuple([p[1]*q[2] - p[2]*q[1],
                             p[2]*q[0] - p[0]*q[2],
                             p[0]*q[1] - p[1]*q[0]])


        modnormal = math.sqrt(self.normal[0]**2 + self.normal[1]**2 + self.normal[2]**2)
        if (modnormal < 1e-15):
            sys.exit("Erro ao calcular o vetor normal")

        self.normal = tuple([self.normal[0]/modnormal, self.normal[1]/modnormal, self.normal[2]/modnormal])

        # salva o parent para pedir a temperatura dos vizinhos
        self.parent = parent

    def link(self, left, right):
        # Estabelece as pastilhas vizinhas
        self.left = left
        self.right = right

    def calc_temp(self):
        # Calcula a teperatura no próximo timestep
        def perimeter_temp(self):
                # obtém a temperatura média dos aneis vizinhos
                bottom_top_ring_temp = self.parent.neighborhood()
                # faz a média ponderada
                return ((self.left.last_temp * self.dl) +
                        (self.right.last_temp * self.dl) +
                        (bottom_top_ring_temp[0] * self.d) +
                        (bottom_top_ring_temp[1] * self.d)) / (self.dl * 2 + self.d * 2)

        self.t += 1
        if self.bursted:
            # se estourou, é só rejunte e a temperatura é a média dos vizinhos
            self.new_temp = perimeter_temp(self)
        else:
            # v . n
            vn = math.fsum([x * y for x, y in zip (self.vel, self.normal)])

            # atrito
            if vn > 0:
                delta_atrito = self.alpha * vn * math.atan( (self.t - self.t_0)**2 )
            else:
                delta_atrito = 0

            # dissipação
            delta_diss = self.delta * abs(vn)

            # nova temperatura
            self.new_temp = perimeter_temp(self) + delta_atrito - delta_diss

            # caso atinja o limite
            if self.new_temp > self.max_temp:
                self.bursted = True
                self.new_temp = (self.left.last_temp + self.right.last_temp) / 2
            

    def update_temp(self):
        # avança com o timestep 
        self.last_temp = self.new_temp
        return self.last_temp

    def __str__(self):
        # formata a impressão
        if self.bursted:
            return "%.2f" % -self.last_temp
        else:
            return "%.2f" % self.last_temp	

class ring():
    # classe dos anéis
    def __init__(self, inputs, l, L):

        def n_tiles(l):
            # Calcula o número de pastilhas em cada anel
            return (2 * math.pi * math.sqrt(l / inputs["a"])) / inputs["d"]

        # lista que terá os aneis
        self.tiles = []

        # Z do circulo inferior (0) e superior(1)
        z0 = l
        z1 = l+L

        # numero de pastilhas (float)
        n = int(n_tiles(l))
        
        # raios
        r0 = math.sqrt(z0 / inputs["a"])
        r1 = math.sqrt(z1 / inputs["a"])

        # angulo de cada pastilha
        alpha0 = 2 * math.asin((inputs["d"]/2) / r0 )
        alpha1 = 2 * math.asin((inputs["d"]/2) / r1 ) 

        # angulo da direfença entre cada pastilha (anel inferior e superior)
        beta0 = ((2 * math.pi) - (n * alpha0)) / n
        beta1 = ((2 * math.pi) - (n * alpha1)) / n 

        for t in range(int(n)):
            # gera cada uma das pastilhas
            x0 = r0 * math.cos( t * (alpha0 + beta0) )
            y0 = r0 * math.sin( t * (alpha0 + beta0) )

            X0 = r0 * math.cos( (t+1) * (alpha0 + beta0) )
            Y0 = r0 * math.sin( (t+1) * (alpha0 + beta0) )

            x1 = r1 * math.cos( t * (alpha1 + beta1) )
            y1 = r1 * math.sin( t * (alpha1 + beta1) )

            X1 = r1 * math.cos( (t+1) * (alpha1 + beta1) )
            Y1 = r1 * math.sin( (t+1) * (alpha1 + beta1) )
            
            self.tiles.append(tile(inputs, (x0, y0, z0),
                                           (X0, Y0, z0),
                                           (x1, y1, z1),
                                           (X1, Y1, z1), self))

        # conecta os vizinhos
        for t in range(len(self.tiles)):
            try:
                self.tiles[t].link(self.tiles[t-1], self.tiles[t+1])
            except IndexError:
                # caso seja a última, liga com o primeiro
                self.tiles[t].link(self.tiles[t-1], self.tiles[0])

        # sua propria temperatura
        self.temp = inputs["theta_0"]
        # referência para os aneis vizinhos
        self.previus_ring = False
        self.next_ring = False
        
    def calc_temp(self):
        # requisita que cada pastilha gere sua próxima temperatura
        for t in self.tiles:
            t.calc_temp()

    def update_temp(self):
        # atualiza a temperatura das pastilhas e sua média
        s = 0
        for t in self.tiles:
            s += t.update_temp()
        # temperatura média
        self.temp = s / len(self.tiles)

    def neighborhood(self):
        # retorna a temperatura dos aneis vizinhos
        return (self.previus_ring.temp, self.next_ring.temp)

    def __str__(self):
        # formata a impressão
        s = ''
        for t in self.tiles:
            s += '%s ' % t
        #s += '\n'
        return s

class cover():
    # calota
    # misto das classes de anel e pastilha, por isso não derivei. Mas acho que deveria
    def __init__(self, inputs):
        self.t_0 = inputs["t_0"]
        self.t = self.t_0

        self.alpha = inputs["alpha"]
        self.delta = inputs["delta"]

        self.last_temp = inputs["theta_0"]
        self.max_temp = inputs["theta_crit"]

        self.vel = inputs["vel"]
        self.pos = inputs["pos"]

        self.bursted = False # True quando a pastilha estourou quando estourou

        # Gabriel: a normal é o inverso do vetor posição (certo?)
        # Max: Estamos considerando que a normal está na mesma direção e sentindo do vetor posição
        #self.normal = tuple(map(lambda x: -x, inputs["pos"]))
        self.normal = inputs["pos"]

        # primeiro anel
        self.next_ring = False
        self.temp = inputs["theta_0"]

    def calc_temp(self):
        self.t += 1
        if self.bursted:
            # se estourou, é só rejunte e a temperatura é a média dos vizinhos
            self.new_temp = self.next_ring.temp
        else:
            # v . n
            vn = math.fsum([x * y for x, y in zip (self.vel, self.normal)])

            # atrito
            if vn > 0:
                delta_atrito = self.alpha * vn * math.atan( (self.t - self.t_0)**2 )
            else:
                delta_atrito = 0

            # dissipação
            delta_diss = self.delta * abs(vn)

            self.new_temp = self.next_ring.temp + delta_atrito - delta_diss

            if self.new_temp > self.max_temp:
                self.bursted = True
                self.new_temp = self.next_ring.temp

    def update_temp(self):
        self.temp = self.last_temp = self.new_temp
        return self.last_temp

    def __str__(self):
        if self.bursted:
            return "%.2f" % -self.last_temp
        else:
            return "%.2f" % self.last_temp	


class mesh():
    # classe que junta todos os aneis mais calotas.
    def __init__(self, inputs):
        self.rings = []

        # L -> altura da menor circunferência do primeiro anel
        L = inputs["a"] * (3 * inputs["d"] / math.pi) ** 2
        l = L

        # calota
        self.cover = cover(inputs)

        # cria os aneis, já referenciando ao anterior
        prev_ring = self.cover
        while ((l + L) < inputs["h"]):
            r = ring(inputs, l, L)
            r.previus_ring = prev_ring
            self.rings.append(r)
            prev_ring = r
            l += L

        # liga os vizinhos "next"
        next_ring = self.rings[-1]
        for r in reversed(self.rings):
            r.next_ring = next_ring
            next_ring = r

        # liga a calota ao primeiro anel
        self.cover.next_ring = self.rings[0]


    def step(self):
        # calcula as novas temperaturas em todos os pontos
        for r in self.rings:
            r.calc_temp()
        self.cover.calc_temp()

        # salva os valores para o próximo step
        for r in self.rings:
            r.update_temp()
        self.cover.update_temp()


def render(mesh, inputs):
    # renderização com o VTK
    def vtkTile(tile, min_temp, max_temp):
        # calcula a intensidade da cor para representar a temperatura.
        # vermelho, com pastilha. azul, estourou
        intensity = 255 - int(((tile.last_temp - min_temp) / (max_temp - min_temp)) * 256)
        if tile.bursted:
            color = (intensity, intensity, 255)
        else:
            color = (255, intensity, intensity)

        Points = vtk.vtkPoints()
        Triangles = vtk.vtkCellArray()

        Points.InsertNextPoint(*tile.edges[0])
        Points.InsertNextPoint(*tile.edges[1])
        Points.InsertNextPoint(*tile.edges[2])
        Points.InsertNextPoint(*tile.edges[3])

        Triangle1 = vtk.vtkTriangle();
        Triangle1.GetPointIds().SetId(0, 0);
        Triangle1.GetPointIds().SetId(1, 1);
        Triangle1.GetPointIds().SetId(2, 2);
        Triangles.InsertNextCell(Triangle1);

        Triangle2 = vtk.vtkTriangle();
        Triangle2.GetPointIds().SetId(0, 1);
        Triangle2.GetPointIds().SetId(1, 3);
        Triangle2.GetPointIds().SetId(2, 2);
        Triangles.InsertNextCell(Triangle2);

        Colors = vtk.vtkUnsignedCharArray();
        Colors.SetNumberOfComponents(3);
        Colors.SetName("Colors");
        Colors.InsertNextTuple3(*color);
        Colors.InsertNextTuple3(*color);

        polydata = vtk.vtkPolyData()
        polydata.SetPoints(Points)
        polydata.SetPolys(Triangles)
 
        polydata.GetCellData().SetScalars(Colors);
        polydata.Modified()
        polydata.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInput(polydata)
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetAmbient(1.0)
        actor.GetProperty().SetDiffuse(0.0)

        return actor

    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
 
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # calcula as temperaturas maximas e minimas para a escala de cor
    # longe do ótimo, mas funciona
    temps = []    

    for ring in mesh.rings:
        for tile in ring.tiles:
            temps.append(tile.last_temp)

    min_temp = min(temps)
    max_temp = max(temps)

    # assign actor to the renderer
    for ring in mesh.rings:
        for tile in ring.tiles:
            ren.AddActor(vtkTile(tile, min_temp, max_temp))

    #ren.SetBackground(1,1,1) # Background color white
 
    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

# rotaciona sistema
def rotate_system(inputs):
    pos = list(inputs["pos"])
    vel = list(inputs["vel"])

    modpos = math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2)
  
    if (modpos < 1e-15):
      sys.exit("Vetor inválido")
  
    pos = [pos[0]/modpos, pos[1]/modpos, pos[2]/modpos]

    if (pos[0] == 0):
      alpha = math.pi / 2.0
    else:
      alpha = (-1) * math.atan( pos[1] / pos[0] )
    
 
    x_linha = vel[0] * math.cos(alpha) - vel[1] * math.sin(alpha)
    y_linha = vel[0] * math.sin(alpha) + vel[1] * math.cos(alpha)
 
    vel[0] = x_linha
    vel[0] = y_linha    

    # rotação da posição
    # rotacao em relacao ao eixo z
    x_linha = pos[0] * math.cos(alpha) - pos[1] * math.sin(alpha)
    y_linha = pos[0] * math.sin(alpha) + pos[1] * math.cos(alpha)
 
    pos[0] = x_linha
    pos[1] = y_linha
 
    # rotacao em relacao ao eixo x (a fim de zerar y)
    if (pos[2] == 0):
       beta = math.pi / 2.0
    else:
       beta = (-1) * math.atan( pos[0] / pos[2])
  
    x_linha = vel[2] * math.sin(beta) + vel[0] * math.cos(beta)
    z_linha = vel[2] * math.cos(beta) - vel[0] * math.sin(beta) 
       
    vel[0] = x_linha
    vel[2] = z_linha
 
    # rotacao em relacao ao eixo x
    x_linha = pos[2] * math.sin(beta) + pos[0] * math.cos(beta)
    z_linha = pos[2] * math.cos(beta) - pos[0] * math.sin(beta) 
       
    pos[0] = x_linha
    pos[2] = z_linha

    if (pos[2] > 0.):
      pos = (-pos[0], -pos[1], -pos[2])
      vel = (-vel[0], -vel[1], -vel[2])
 
    inputs["pos"] = (pos[0], pos[1], pos[2])
    inputs["vel"] = (vel[0], vel[1], vel[2])

    # Estamos considerando que o vetor posição está saindo da cápsula e não entrando


def parse_input():
    i = {}
    
    r = re.compile('[ \t\n\r]+')
    values = map(float, r.split(sys.stdin.read())[:-1])
    i["h"], i["a"] = values[0:2]
    i["d"] = values[2]
    i["alpha"], i["t_0"] = values[3:5]
    i["delta"] = values[5]
    i["theta_crit"], i["theta_0"] = values[6:8]
    i["pos"] = tuple(values[8:11])
    i["vel"] = tuple(values[11:14])
    i["steps"] = int(values[14])
    
#    i["h"] = 450
#    i["a"] = 3
#    i["d"] = 3.5
#    i["alpha"] = 2
#    i["delta"] = 0.5
#    i["t_0"] = 0.0
#    i["t_inicial"] = 0
#    i["theta_crit"] = 1000 #?
#    i["theta_0"] = 10
#    i["pos"] = (1.0, 1.0, 1.0)
#    i["vel"] = (0, 0, -5)
#    i["steps"] = 1000
    return i

if __name__ == "__main__":
    # le o arquivo de entrada (ou deveria...)
    inputs = parse_input()
    # gera a malha
    m = mesh(inputs)
    # executa todos os steps
    for ts in range(inputs["steps"]):
        m.step()

    # gera a saida
    print "%.2f %.2f %.2f" % (inputs["a"], inputs["h"], inputs["d"])
    print m.cover
    print "%d" % len(m.rings)
    for ring in m.rings:
        print ring

    # exibe em 3D
    if vtk:
        render(m, inputs)

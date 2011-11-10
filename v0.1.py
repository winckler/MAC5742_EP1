#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Esboco do EP de MAC5742 - 2/2011
# 
# Copyright (C) 2011 Gabriel A. von Winckler
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

# version 0.1
#
# Bugs:
#  - Cálculo MUITO simplificado, ainda não aprovado pelo Gubi !!!
#  - A normal está calculada nas coordenadas do elipsoide

import math

class tile():
    def __init__(self, inputs, a, b, c, d):

        self.t_0 = inputs["t_0"]
        self.t = self.t_0

        self.alpha = inputs["alpha"]
        self.delta = inputs["delta"]

        self.last_temp = inputs["theta_0"]
        self.max_temp = inputs["theta_crit"]

        self.vel = inputs["vel"]
        self.pos = inputs["pos"]

        self.bursted = False # True quando a pastilha estourou quando estourou

        # produto vetorial n = (p2 - p1) x (p3 - p1)
        p = [x - y for x, y in zip (b, a)]
        q = [x - y for x, y in zip (d, a)]

        # normal na coordenada do paraboloide. CORRIGIR!
        self.normal = tuple([p[1]*q[2] - p[2]*q[1],
                             p[2]*q[0] - p[0]*q[2],
                             p[0]*q[1] - p[1]*q[0]])

    def link(self, left, right):
        self.left = left
        self.right = right

    def calc_temp(self):
        self.t += 1
        if self.bursted:
            # se estourou, é só rejunte e a temperatura é a média dos vizinhos
            self.new_temp = (self.left.last_temp + self.right.last_temp) / 2
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

            self.new_temp = ((self.left.last_temp + self.right.last_temp) / 2) + delta_atrito - delta_diss

            if self.new_temp > self.max_temp:
                self.bursted = True
                self.new_temp = (self.left.last_temp + self.right.last_temp) / 2
            

    def update_temp(self):
        self.last_temp = self.new_temp

    def __repr__(self):
        if self.bursted:
            return "%.2f" % -self.last_temp
        else:
            return "%.2f" % self.last_temp	

def parse_input():
    # Não está lendo o do arquivo!
    i = {}
    i["h"] = 22.0
    i["a"] = 2.0
    i["d"] = 2.0
    i["alpha"] = 1.0
    i["delta"] = 1.0
    i["t_0"] = 0.0
    i["t_inicial"] = 0
    i["theta_crit"] = 2000 #?
    i["theta_0"] = 1000.0 # really?
    i["pos"] = (1.0, 1.0, 1.0)
    i["vel"] = (1.0, 1.0, 1.0)
    i["steps"] = 10
    return i

def generate_mesh(i):
    mesh = []

    def n_tiles(l):
        return (2 * math.pi * math.sqrt(l / i["a"])) / i["d"]

    # L -> altura da menor circunferência do primeiro anel
    L = i["a"] * (3 * i["d"] / math.pi) ** 2
    l = L
    # Cada anel
    print l+L
    while ((l + L) < i["h"]):
        ring = []

        # Z do circulo inferior (0) e superior(1)
        z0 = l
        z1 = l+L

        # numero de pastilhas (float)
        n = int(n_tiles(l))
        
        # raios
        r0 = math.sqrt(z0 / i["a"])
        r1 = math.sqrt(z1 / i["a"])

        # angulo de cada pastilha
        alpha0 = 2 * math.asin((i["d"]/2) / r0 )
        alpha1 = 2 * math.asin((i["d"]/2) / r1 ) 

        # angulo da direfença entre cada pastilha (anel inferior e superior)
        beta0 = ((2 * math.pi) - (n * alpha0)) / n
        beta1 = ((2 * math.pi) - (n * alpha1)) / n 

        for t in range(int(n)):
            
            x0 = r0 * math.cos( t * (alpha0 + beta0) )
            y0 = r0 * math.sin( t * (alpha0 + beta0) )

            X0 = r0 * math.cos( (t+1) * (alpha0 + beta0) )
            Y0 = r0 * math.sin( (t+1) * (alpha0 + beta0) )

            x1 = r1 * math.cos( t * (alpha1 + beta1) )
            y1 = r1 * math.sin( t * (alpha1 + beta1) )

            X1 = r1 * math.cos( (t+1) * (alpha1 + beta1) )
            Y1 = r1 * math.sin( (t+1) * (alpha1 + beta1) )
            
            ring.append(tile(inputs,(x0, y0, z0),
                                    (X0, Y0, z0),
                                    (x1, y1, z1),
                                    (X0, Y0, z1))) 

        # conectar os vizinhos
        for t in range(len(ring)):
            try:
                ring[t].link(ring[t-1], ring[t+1])
            except IndexError:
                ring[t].link(ring[t-1], ring[0])

        mesh.append(ring)
        l += L

    return mesh

def run(mesh, inputs):
    # calcula as novas temperaturas em todos os pontos
    for ring in mesh:
        for tile in ring:
            tile.calc_temp()

    # salva os valores para o próximo step
    for ring in mesh:
        for tile in ring:
            tile.update_temp()

def output(mesh, inputs):
    print "%.2f %.2f %.2f" % (inputs["a"], inputs["h"], inputs["d"])
    print "%.2f (nao calculado)" % (inputs["theta_0"])
    print "%d" % len(mesh)
    for ring in mesh:
        for tile in ring:
            print tile,
        print

if __name__ == "__main__":
    inputs = parse_input()
    mesh = generate_mesh(inputs)
    for ts in range(inputs["steps"]):
        run(mesh, inputs) 
    output(mesh, inputs)

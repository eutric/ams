function [O_best, th_best_12] = tappeto(O1, O2, fun_ob, n, a, b)
% Funzione che ottimizza la function handle data come input, tra a e b
% provando ogni combinazione possibile del dominio discretizzato.
% INPUT
% O1
% OUTPUT
%
costmin = 1e25;

tic
for th1i = linspace(a(1), b(1), n)
    for th2f = linspace(a(2), b(2), n)
        for om = linspace(a(3), b(3), n)
            x_t = [th1i, th2f, om];
            costo = fun_ob(x_t);
            O_t = O_tfun(O1, O2, x_t(1), x_t(2), x_t(3), [0,0]);
            if O_t.e<1 && O_t.e>=0
                if costo < costmin
                    O_best = O_t;
                    O_best.x = x_t;
                    th_best_12 = [th1i, th2f];
                end
            end
        end
    end
end
toc

% tic
% for th1i=linspace(a, b, n)
%     for th2f=linspace(a, b, n)
%         for om=linspace(1, 2*pi,n)
%             [O_t] = O_tfun(O_start,O_end,th1i,th2f,om, [0,0]);
%             if O_t.e<1 && O_t.e>=0
%                 OO_t=[OO_t,O_t];
%                 if O_t.cost<costmin
%                     costmin=O_t.cost;
%                     O_best=O_t;
%                     th_t_best=O_t.th_t;
%                     th_best_12 = [th1i, th2f];
%                 end
%             else 
%                 kk = kk+1;
%             end
% 
%         end
%     end
% end
% toc

end


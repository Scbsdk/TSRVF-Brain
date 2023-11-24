function BestK = estimate(D)

      Sim = exp(-(D/4).^2);
      S_eps = Sim; S_eps(S_eps<0.5) = 0;
      G_eps = graph(S_eps);
      for i=1:10
          BestK = i;
          if(i ==  max(unique(conncomp(G_eps)))) break; 
          end
      end
      
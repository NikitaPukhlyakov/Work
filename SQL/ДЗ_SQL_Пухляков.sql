@set date_now =  current_date()
@set date_start = current_date()  - interval 1 year
CREATE TABLE analytics_temp.pukhlyakov
(  PRIMARY KEY (id),
   order_id int,
   product varchar(20),
   tickets_count int,
   email varchar(255),
   phone varchar(127),
   profit decimal (9,2),
   device_type varchar(20),
   site_version varchar(30),
   os varchar(300), 
   cdate datetime
);
-- CREATE TABLE analytics_temp.pukhlyakov_ryobject
-- (  id_key int NOT NULL AUTO_INCREMENT,
--    id int,
--    arrival_station_id int,
--    PRIMARY key(id_key)
-- );

-- avia
INSERT INTO analytics_temp.pukhlyakov(order_id,
                                      product,
                                      tickets_count,
   									  email, 
                                      phone, 
                                      profit, 
                                      device_type, 
                                      site_version,
                                      os,
                                      cdate)
SELECT B.id, 'avia', B.tickets_count, B.contact_email,B.contact_phone, S.tutu_order_revenue, I.device_type , I.site_version ,I.os , B.cdate
FROM favia.rm_forder_object_backoffice B
INNER JOIN favia.rm_forder_object_sales S ON B.id = S.id 
INNER JOIN favia.rm_forder_object_interface_statistic I ON B.id=I.id
WHERE B.is_sold = 1
           AND B.contact_email NOT LIKE '%@tutu.ru%'
           AND B.gds != 'fake'
           AND (B.cdate BETWEEN :date_start AND :date_now);
-- bus
INSERT INTO analytics_temp.pukhlyakov (order_id,
                                       product,
                                       tickets_count,
                                       email,
                                       phone,
                                       device_type,
                                       site_version,
                                       profit,
                                       os,
                                       cdate)
SELECT BO.id, 'bus', BO.passengers_count, BO.email ,BO.phone, BA.device_type,BA.site_version, (BT.price - BT.price_original) * 0.5 * (1-18/118) ,BA.os , BO.cdate
FROM busorder.order_backoffice BO
INNER JOIN busorder.analytics_data BA  ON BO.id = BA.order_id
INNER JOIN busorder.ticket BT ON  BO.id = BT.order_id 
WHERE BO.email NOT LIKE '%@tutu.ru'
           AND BO.state IN ('issued_authorized', 'issued_not_authorized')
           AND ( BO.cdate BETWEEN :date_start AND :date_now );
 --    rzd
 INSERT INTO analytics_temp.pukhlyakov (order_id,
                                        product,
                                        tickets_count,
                                        email,
                                        phone,
                                        #device_type,
                                        #site_version,
                                        profit,
                                        #os
                                        cdate)
SELECT RYO.id, 'rzd', RYO.paid_tickets, RYO.email,RYO.phone,RYO.profit, RYO.cdate
FROM avia.rm_ryobject RYO
WHERE RYO.state = 'paid' 
           AND processed_by NOT IN ('fake', 'fakeUit')
           AND (RYO.cdate BETWEEN :date_start AND :date_now );
          
     
-- tour
INSERT INTO analytics_temp.pukhlyakov (order_id,
                                       product,
                                       tickets_count,
                                       email,
                                       phone,
                                       device_type,
                                       site_version,
                                       profit,
                                       os,
                                       cdate)
SELECT TOR.id,'tour', TOR.tickets_included, OC.email,OC.phone, NULL, NULL, (PR.customer_total_price_rub - ACC.net_cost_rub) * (1-18/118) - 0.0125 * PR.customer_total_price_rub , NULL, TOR.cdate
FROM avia.rm_tour_order TOR
INNER JOIN avia.rm_tour_order_price PR ON TOR.id=PR.order_id
INNER JOIN avia.rm_tour_accounting_data ACC ON TOR.id=ACC.id
INNER JOIN avia.rm_tour_order_customer OC ON TOR.id=OC.order_id
        WHERE status IN ('confirmed') AND is_external=0
        AND (TOR.cdate BETWEEN :date_start AND :date_now )
        AND NOT OC.email IN  ('tolstaya@tutu.ru','tours@tutu.ru','doctor_wh0@inbox.ru', 'test@test.com', 'rodionova@tutu.ru', 'pinkas@tutu.ru', 'lipatov@tutu.ru', 'nbelan@tutu.ru', 'qatest@tutu.ru');                            
       
-- подгрузка жд       
UPDATE analytics_temp.pukhlyakov, avia.rm_ryorder_analytics_data 
SET analytics_temp.pukhlyakov.device_type = avia.rm_ryorder_analytics_data.device_type
WHERE analytics_temp.pukhlyakov.order_id = avia.rm_ryorder_analytics_data.object_id
      AND analytics_temp.pukhlyakov.product='rzd'; 
     
UPDATE analytics_temp.pukhlyakov, avia.rm_ryorder_analytics_data 
SET analytics_temp.pukhlyakov.os =  (case when avia.rm_ryorder_analytics_data.device_type!='PC' then SUBSTRING_INDEX(avia.rm_ryorder_analytics_data.user_from,'-',-1)  else NULL end )      
WHERE analytics_temp.pukhlyakov.order_id = avia.rm_ryorder_analytics_data.object_id
      AND analytics_temp.pukhlyakov.product='rzd'; 
-- добавил site_version
UPDATE analytics_temp.pukhlyakov, avia.rm_ryorder_analytics_data 
SET analytics_temp.pukhlyakov.site_version = (case when avia.rm_ryorder_analytics_data.browser 
                                                        NOT IN ('app','mobile','tuturu','NULL')
                                                   then 'browser'
                                                   else  avia.rm_ryorder_analytics_data.browser
                                                   end)
WHERE analytics_temp.pukhlyakov.order_id = avia.rm_ryorder_analytics_data.object_id
      AND analytics_temp.pukhlyakov.product='rzd';     
-- телефоны      
UPDATE analytics_temp.pukhlyakov 
SET analytics_temp.pukhlyakov.phone = (case when analytics_temp.pukhlyakov.phone like '8%'
                                             then '+7' + SUBSTRING(analytics_temp.pukhlyakov.phone,2,10)
                                             else  analytics_temp.pukhlyakov.phone
                                        end)
WHERE LENGTH(analytics_temp.pukhlyakov.phone)=11;   

-- email
UPDATE analytics_temp.pukhlyakov 
SET analytics_temp.pukhlyakov.email = LOWER(analytics_temp.pukhlyakov.email);
     
-- QUESTIONS
CREATE INDEX Clientgroup ON analytics_temp.pukhlyakov (phone); 
-- так как группируем по полю phone

-- a
-- Каков % людей (и какой % профита они нам принесли), у которых есть покупки во всех 4х продуктах? (Использовать having + group by + distinct)

SELECT COUNT(A.phone)/(SELECT COUNT(DISTINCT phone) FROM analytics_temp.pukhlyakov)*100, SUM(A.profit)/(SELECT SUM(profit) FROM analytics_temp.pukhlyakov)*100
FROM(
     SELECT phone, SUM(profit) as profit
	 FROM analytics_temp.pukhlyakov
	 GROUP BY phone
	 HAVING COUNT(DISTINCT product)=4
	 ) A;

-- b
-- Аналогично для 3х продуктов (выберите любые) (Нужно расставить в таблице индексы, сделать join и можно использовать дополнительные таблицы или temporary table)	
SELECT COUNT(A.phone)/(SELECT COUNT(DISTINCT phone) FROM analytics_temp.pukhlyakov)*100, SUM(A.profit)/(SELECT SUM(profit) FROM analytics_temp.pukhlyakov)*100
FROM(
     SELECT phone, SUM(profit) as profit
	 FROM analytics_temp.pukhlyakov
	 WHERE analytics_temp.pukhlyakov.product IN ('avia','rzd','bus')
	 GROUP BY phone
	 HAVING COUNT(DISTINCT product)=3
	 ) A;

-- c
-- Продукт, приносящий наибольший профит (и его % от суммарного профита)? (Нужно быть готовым посчитать "на лету" профит для других продуктов, необходимо использовать group by)	
SELECT product, SUM(profit)/(SELECT SUM(profit) FROM analytics_temp.pukhlyakov)*100 AS percent
FROM analytics_temp.pukhlyakov
GROUP BY product
ORDER BY 2 DESC;

-- d
-- % заказов в ЖД, у людей, у которых нет авиа заказов? 
SELECT COUNT(A.phone)/(SELECT COUNT(phone) FROM analytics_temp.pukhlyakov WHERE product='rzd')*100 AS percent
FROM(SELECT phone  
     FROM analytics_temp.pukhlyakov
	 WHERE product='rzd'
	 GROUP BY phone
	 HAVING MIN(product)>'avia'
	 ) A;
 -- не повис:)

-- e
-- Определить людей, сделавших наибольшее кол-во покупок в каждом продукте и результат вывести в одну таблицу и отсортировать по количеству покупок	
(SELECT phone, COUNT(phone) as num, product
FROM analytics_temp.pukhlyakov
WHERE product='avia'
GROUP BY phone
ORDER BY 2 DESC
LIMIT 1)
UNION
(SELECT phone, COUNT(phone) as num , product
FROM analytics_temp.pukhlyakov
WHERE product='bus'
GROUP BY phone
ORDER BY 2 DESC
LIMIT 1)
UNION
(SELECT phone, COUNT(phone) as num, product
FROM analytics_temp.pukhlyakov
WHERE product='rzd'
GROUP BY phone
ORDER BY 2 DESC
LIMIT 1)
UNION
(SELECT phone, COUNT(phone) as num , product
FROM analytics_temp.pukhlyakov
WHERE product='tour'
GROUP BY phone
ORDER BY 2 DESC
LIMIT 1)
ORDER BY 2 DESC;

-- f
-- % людей, делавших первую покупку в Автобусах, а потом, в течение года, делавших покупки в других продуктах? 
SELECT COUNT(DISTINCT C.phone)/(SELECT COUNT(DISTINCT phone) FROM analytics_temp.pukhlyakov)*100 AS percent
FROM
(SELECT A.phone as phone, P.product as product
FROM(
	   SELECT phone, MIN(cdate) as min_date
	   FROM analytics_temp.pukhlyakov
	   GROUP BY phone
	   HAVING COUNT(DISTINCT product)>1
	)A   
	   INNER JOIN  analytics_temp.pukhlyakov P
	              ON (A.min_date=P.cdate and A.phone=P.phone)
)C
WHERE C.product='bus'; 
-- g
-- Сколько людей пользуется почтой Яндекс? (и в %) 
SELECT COUNT(DISTINCT phone), COUNT(DISTINCT phone)/(SELECT COUNT(DISTINCT phone) FROM analytics_temp.pukhlyakov)*100 AS percent
FROM analytics_temp.pukhlyakov
WHERE email LIKE '%@yandex.ru'

-- h
-- По месяцам по всем продуктам (по турам не смотрим, у туров нет этой информации) посчитать для каждой из версий сайта (и потом для каждой из ОС), сколько людей делали заказы
-- для версии сайта
SELECT product, MONTHNAME(cdate) as month_ ,
           count(distinct case when site_version='mobile' then phone else null  end) 'mobile',
           count(distinct case when site_version='app' then phone else null  end) 'app',
           count(distinct case when site_version='browser' then phone else null  end) 'browser',
           count(distinct case when site_version='tuturu' then phone else null  end) 'tuturu'
FROM analytics_temp.pukhlyakov
WHERE product NOT IN ('tour')
GROUP BY product, MONTHNAME(cdate);


-- для ОС
SELECT product, MONTHNAME(cdate) as month_ , 
           count(distinct case when os LIKE '%ios%' then phone else null  end) 'ios',
           count(distinct case when os LIKE '%ndroid%' then phone else null  end) 'android',
           count(distinct case when os LIKE '%indows%'then phone else null  end) 'windows'
FROM analytics_temp.pukhlyakov
WHERE product NOT IN ('tour')
	      AND device_type!='PC'
GROUP BY product, MONTHNAME(cdate);

